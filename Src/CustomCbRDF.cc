/*
 *  Adyton: A Network Simulator for Opportunistic Networks
 *  Copyright (C) 2015-2018  Nikolaos Papanikos,
 *  Dimitrios-Georgios Akestoridis, and Evangelos Papapetrou
 *
 *  This file is part of Adyton.
 *
 *  Adyton is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Adyton is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Adyton.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Written by Evangelos Papapetrou and Nikolaos Papanikos.
 */


#ifndef CUSTOMCustomCbRDF_H
#define CUSTOMCustomCbRDF_H

#include "CustomCbRDF.h"

#endif

//#define CustomCbRDF_DEBUG

// +----------------+
// | Protocol steps |
// +----------------+----------------------------------------------------------+
// | Node (a) encounters node (b)                                              |
// |===========================================================================|
// | (a): [method: Contact()]                                                  |
// |      - Updates the value of the utility in use.                           |
// |===========================================================================|
// | <steps related to direct delivery are performed (check Routing.cc)>       |
// |===========================================================================|
// | (a): [method: AfterDirectTransfers()]                                     |
// |      - Checks the utility given by the user. If PRoPHET or SimBetTS       |
// |        is selected then the optional steps that form the ego/contact      |
// |        graph are executed.                                                |
// |      +--------------------------------------------+                       |
// |      | Formation of ego/contacts graph (optional) |                       |
// |      +--------------------------------------------+-----------------+     |
// |      | (a): [method: SendContactRequest()]                          |     |
// |      |      - Sends an empty packet to request the direct contacts  |     |
// |      |        of node (b).                                          |     |
// |      |                                                              |     |
// |      | (a)--------------contact request-------------->(b)           |     |
// |      |                                                              |     |
// |      | (b): [method: ReceptionRequestContacts()]                    |     |
// |      |      - Creates a packet that contains the node's direct      |     |
// |      |        contacts. In case of PRoPHET an array with the        |     |
// |      |        delivery probability of node (b) to reach all         |     |
// |      |        destinations is piggybacked to the packet.            |     |
// |      |                                                              |     |
// |      | (a)<-------------contact/DPT reply-------------(b)           |     |
// |      |                                                              |     |
// |      | (a): [method: ReceptionContacts()] - for contact reply -     |     |
// |      |      - Updates its ego/contact data structure using the      |     |
// |      |        <contact reply>.                                      |     |
// |      | (a): [method: ReceptionDPT()] - for DPT reply -              |     |
// |      |      - Updates its DPT data structure using the              |     |
// |      |        <DPT reply>.                                          |     |
// |      +--------------------------------------------------------------+     |
// |      (a): [method: (SendSummary())]                                       |
// |      - Creates a <summary packet> that contains a list of packet IDs that |
// |        node (a) has for node (b). Also, the destination ID of each packet |
// |        is piggybacked to this <summary packet>.                           |
// |                                                                           |
// |      (a)----summary packet (packet IDs, Destinations)---->(b)             |
// |                                                                           |
// |===========================================================================|
// | (b): [method: ReceptionSummary()]                                         |
// |      - For each packet ID inside the <summary packet> node (b):           |
// |            * Checks if the packet exists in its buffer.                   |
// |            * Calculates the utility value that it has for the packet.     |
// |      - Creates a <request packet> that contains the packet IDs along      |
// |        with the utility value that node (b) has for each packet.          |
// |                                                                           |
// | (a)<------------request packet (packet IDs, Utils)-----------------(b)    |
// |                                                                           |
// |===========================================================================|
// | (a): [method: ReceptionRequest()]                                         |
// |      - Gets the requested packets from its buffer.                        |
// |      - For each packet that node (b) does not have:                       |
// |            * If the node's (b) utility value is larger than "tau" create  |
// |              a new packet replica and transmit.                           |
// |                                                                           |
// |     ** "tau": largest utility value of a network node that node (a) has   |
// |               transmitted a packet replica to.                            |
// |                                                                           |
// | (a)------------------------packet replica 1 ---------------------->(b)    |
// | (a)------------------------packet replica 2 ---------------------->(b)    |
// |                                  ...                                      |
// | (a)------------------------packet replica n ---------------------->(b)    |
// |                                                                           |
// |===========================================================================|
// |(b): [method: ReceptionData()]                                             |
// |     - Receives the transmitted packets one by one.                        |
// |     - Sets the "tau" value equal to node's (b) utility value.             |
// |===========================================================================|
// | (a): [method: ContactRemoved()]                                           |
// |      - Updates the value of the utility in use.                           |
// +---------------------------------------------------------------------------+


CustomCbRDF::CustomCbRDF(PacketPool *PP, MAC *mc, PacketBuffer *Bf, int NID, Statistics *St, Settings *S, God *G)
        : Routing(PP, mc, Bf, NID, St, S, G) {
    string UtilityType;
    double DensityVal;
    double Pinit;
    double Pmin;
    double Pmax;
    double beta;
    double gamma;
    double delta;
    double agingTimeUnit;
    string profileAttribute;
    Util = NULL;
    Adja = NULL;
    MyDPT = NULL;
    //----------------Added by epap----------------------------//
    TwoDOn = false;
    isLVQon = false;
    isKmeansweighted = false;
    isKmeansnormalized = false;
    KmeansUPeriod = 0;
    //---------------------------------------------------------//
    //----------------Added by Vag-----------------------------//
    DstDepend = true;
    //---------------------------------------------------------//
    if (S->ProfileExists() && (UtilityType = S->GetProfileAttribute("Utility")) != "none") {
        //remove all whitespaces
        UtilityType.erase(std::remove_if(UtilityType.begin(), UtilityType.end(), ::isspace), UtilityType.end());
        if (UtilityType == "LTS") {
            Util = new LTS(NID, Set->getNN());
        } else if (UtilityType == "LastContact") {
            //----------------Added by Vag-----------------------------//
            DstDepend = false;
            //---------------------------------------------------------//
            Util = new LastContact(NID, Set->getNN());
        } else if (UtilityType == "DestEnc") {
            //Util=new DestEnc(NID,Set->getNN(),0.85,3600.0*6.0);
            Util = new DestEnc(NID, Set->getNN());
        } else if (UtilityType == "Enc") {
            //----------------Added by Vag-----------------------------//
            DstDepend = false;
            //---------------------------------------------------------//
            //Util=new Enc(NID,Set->getNN(),0.85,3600.0*8.0);
            Util = new Enc(NID, Set->getNN());
        } else if (UtilityType == "AMT") {
            Util = new AMT(NID, Set->getNN());
        } else if (UtilityType == "AIT") {
            Util = new AIT(NID, Set->getNN());
        } else if (UtilityType == "SPM") {
            Util = new SPM(NID, Set->getNN());
        } else if (UtilityType == "SimBet" || UtilityType == "SimBetTS" || UtilityType == "Sim" ||
                   UtilityType == "Bet") {
            //----------------Added by Vag-----------------------------//
            if (UtilityType == "Bet") {
                DstDepend = false;
            }
            //---------------------------------------------------------//
            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("Density")) != "none") {
                DensityVal = atof(profileAttribute.c_str());

                if (DensityVal < 0.0 || DensityVal > 1.0) {
                    printf("Error: Density value must be between 0 and 1. You gave %f\nExiting\n", DensityVal);
                    exit(1);
                }

                //Create contact graph powered with contact aggregation
                if ((profileAttribute = S->GetProfileAttribute("AggregationType")) == "MR") {
                    //printf("Contact graph powered with MR contact aggregation\n");
                    Adja = new Adjacency(NID, Set->getNN(), 2, DensityVal);
                } else {
                    if (profileAttribute == "MF") {
                        Adja = new Adjacency(NID, Set->getNN(), 1, DensityVal);
                    } else {
                        printf("Error: Invalid AggregationType setting\nExiting...");
                        exit(1);
                    }
                }
            } else {
                Adja = new Adjacency(NID, Set->getNN());
            }
        } else if (UtilityType == "Prophet") {
            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT")) != "none") {
                if (profileAttribute == "v1") {
                    MyDPT = new DPTv1(NID, Set->getNN());
                } else if (profileAttribute == "v1.5") {
                    MyDPT = new DPTv1point5(NID, Set->getNN());
                } else if (profileAttribute == "v2") {
                    MyDPT = new DPTv2(NID, Set->getNN());
                } else if (profileAttribute == "v3") {
                    MyDPT = new DPTv3(NID, Set->getNN());
                } else {
                    printf("Error: Unknown version of DPT (%s)\nExiting...", profileAttribute.c_str());
                    exit(1);
                }
            } else {
                MyDPT = new DPTv3(NID, Set->getNN());
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_Pinit")) != "none") {
                Pinit = atof(profileAttribute.c_str());

                if (Pinit < 0.0 || Pinit > 1.0) {
                    printf("Error: Invalid Pinit value (%f)\nExiting...", Pinit);
                    exit(1);
                }

                MyDPT->setPinit(Pinit);
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_Pmin")) != "none") {
                Pmin = atof(profileAttribute.c_str());

                if (Pmin < 0.0 || Pmin > 1.0) {
                    printf("Error: Invalid Pmin value (%f)\nExiting...", Pmin);
                    exit(1);
                }

                MyDPT->setPmin(Pmin);
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_Pmax")) != "none") {
                Pmax = atof(profileAttribute.c_str());

                if (Pmax < 0.0 || Pmax > 1.0) {
                    printf("Error: Invalid Pmax value (%f)\nExiting...", Pmax);
                    exit(1);
                }

                MyDPT->setPmax(Pmax);
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_beta")) != "none") {
                beta = atof(profileAttribute.c_str());

                if (beta < 0.0 || beta > 1.0) {
                    printf("Error: Invalid beta value (%f)\nExiting...", beta);
                    exit(1);
                }

                MyDPT->setBeta(beta);
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_gamma")) != "none") {
                gamma = atof(profileAttribute.c_str());

                if (gamma < 0.0 || gamma > 1.0) {
                    printf("Error: Invalid gamma value (%f)\nExiting...", gamma);
                    exit(1);
                }

                MyDPT->setGamma(gamma);
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_delta")) != "none") {
                delta = atof(profileAttribute.c_str());

                if (delta < 0.0 || delta > 1.0) {
                    printf("Error: Invalid delta value (%f)\nExiting...", delta);
                    exit(1);
                }

                MyDPT->setDelta(delta);
            }

            if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("DPT_agingTimeUnit")) != "none") {
                agingTimeUnit = atof(profileAttribute.c_str());

                if (agingTimeUnit <= 0.0) {
                    printf("Error: Invalid aging time unit value (%f)\nExiting...", agingTimeUnit);
                    exit(1);
                }

                MyDPT->setAgingTimeUnit(agingTimeUnit);
            }
        } else {
            printf("Please select a proper utility!\n");
            exit(1);
        }
    } else {
        //set as default the LTS Utility
        Util = new LTS(NID, Set->getNN());
    }
    this->UType = UtilityType;
    UUpdate = false;
    if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("update")) != "none") {
        if (profileAttribute == "on") {
            UUpdate = true;
        } else {
            if (profileAttribute == "off") {
                UUpdate = false;
            } else {
                printf("Error: Invalid update setting\nExiting...");
                exit(1);
            }
        }
    }
    if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("2D")) != "none") {//Added by epap
        if (profileAttribute == "on") {
            TwoDOn = true;
        } else {
            if (profileAttribute == "off") {
                TwoDOn = false;
            } else {
                printf("Error: Invalid Kmeans 2D setting\nExiting...");
                exit(1);
            }
        }
    }
    if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("LVQ")) != "none") {//Added by epap
        if (profileAttribute == "on") {
            isLVQon = true;
        } else {
            if (profileAttribute == "off") {
                isLVQon = false;
            } else {
                printf("Error: Invalid Kmeans LVQ setting\nExiting...");
                exit(1);
            }
        }
    }
    if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("Kmeans-norm")) != "none") {//Added by epap
        if (profileAttribute == "on") {
            isKmeansnormalized = true;
        } else {
            if (profileAttribute == "off") {
                isKmeansnormalized = false;
            } else {
                printf("Error: Invalid Kmeans normalization setting\nExiting...");
                exit(1);
            }
        }
    }
    if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("Kmeans-period")) != "none") {//Added by epap

        KmeansUPeriod = atoi(profileAttribute.c_str());
        if (KmeansUPeriod < 0) {
            printf("Error: Invalid Kmeans update period setting\nExiting...");
            exit(1);
        }
    }
    if (S->ProfileExists() && (profileAttribute = S->GetProfileAttribute("Kmeans-weighted")) != "none") {//Added by epap
        if (profileAttribute == "on") {
            isKmeansweighted = true;
        } else {
            if (profileAttribute == "off") {
                isKmeansweighted = false;
            } else {
                printf("Error: Invalid Kmeans weighted setting\nExiting...");
                exit(1);
            }
        }
    }
    return;
}


CustomCbRDF::~CustomCbRDF() {
    delete Util;
    delete Adja;
    delete MyDPT;
    return;
}


void CustomCbRDF::NewContact(double CTime, int NID) {
    Contact(CTime, NID);
    return;
}


void CustomCbRDF::Contact(double CTime, int NID) {
// 	printf("%d -> %d\n",this->NodeID,NID);
    if (Adja != NULL) {
        //Update Ego network
        Adja->SetConnection(this->NodeID, NID, CTime);
        //Update data structures that collect data to estimate metrics
        Adja->ContactStart(NID, CTime);
    }
    if (Util != NULL) {
        Util->ContactUp(NID, CTime);
    }
    if (MyDPT != NULL) {
        MyDPT->ContactUp(NID, CTime);
    }
    //Get information about known delivered packets
    int *Information = DM->GetInfo();
    if (Information != NULL) {
        //Create a vaccine information packet
        SendVaccine(CTime, NID, Information);
    } else {
        //Clean buffer using Deletion method (Delivered pkts)
        DM->CleanBuffer(this->Buf);
        if (DM->ExchangeDirectSummary()) {
            SendDirectSummary(CTime, NID);
        } else {
            SendDirectPackets(CTime, NID);
        }
    }
    return;
}


void CustomCbRDF::ContactRemoved(double CTime, int NID) {
    //printf("%d:Contact %d removed\n",this->NodeID,NID);
    if (Adja != NULL) {
        Adja->ContactEnd(NID, CTime);
    }
    if (Util != NULL) {
        Util->ContactDown(NID, CTime);
    }
    if (MyDPT != NULL) {
        MyDPT->ContactDown(NID, CTime);
    }

    if (UUpdate) {
        //Update max util for stored packets
        //printf("%d:Contact with node %d ended at %f\n",this->NodeID,NID,CTime);
        double util = 0.0;
        if (MyDPT && this->UType == "Prophet") {
            util = MyDPT->getDPto(NID, CTime);
        } else if (Adja) {
            if (this->UType == "Bet") {
                util = Adja->getBet();;
            } else if (this->UType == "Sim") {
                util = Adja->getSim(NID);
            } else {
                util = 0.0;
            }
        } else {
            util = Util->get(NID, CTime);;
        }

        if (UType == "Bet" || UType == "Enc") {
            Buf->UpdateThresholdDI(util);
        } else if (UType == "Sim") {
            int NN = this->Set->getNN();
            double *Sims = (double *) malloc(sizeof(double) * NN);
            for (int i = 0; i < NN; i++) {
                Sims[i] = Adja->getSim(i);
            }
            Buf->UpdateThresholdDD(Sims);
            free(Sims);
        } else if (UType == "SimBetTS" || UType == "SimBet") {
            int NN = this->Set->getNN();
            struct SimBetTSmetrics *localMetrics = (struct SimBetTSmetrics *) malloc(
                    sizeof(struct SimBetTSmetrics) * NN);
            for (int i = 0; i < NN; i++) {
                localMetrics[i].Similarity = Adja->getSim(i);
                localMetrics[i].Betweenness = Adja->getBet();
                localMetrics[i].Frequency = Adja->getFreq(i);
                localMetrics[i].Intimacy = Adja->getIntimacy(i);
                localMetrics[i].Recency = Adja->getRecency(i, CTime);
            }
            if (UType == "SimBetTS") {
                Buf->UpdateThresholdSimBetTS(localMetrics);
            } else {
                if (TwoDOn) {//Added by epap
                    Buf->UpdateThresholdSimBet2D(localMetrics);
                } else {
                    Buf->UpdateThresholdSimBet(localMetrics);
                }
            }
            free(localMetrics);
        } else {
            Buf->UpdateThresholdDD(NID, util);
        }
    }
    return;
}


void CustomCbRDF::AfterDirectTransfers(double CTime, int NID) {
    if (MyDPT || Adja) {
        SendContactRequest(CTime, NID);
    } else {
// 		if(UType == "SPM")
// 		{
// 			SendRSPMRequest(CTime,NID);
// 		}
// 		else
// 		{
// 			SendSummary(CTime,NID);
// 		}
        SendSummary(CTime, NID);
    }
    return;
}


void CustomCbRDF::recv(double rTime, int pktID) {
    int PacketID = 0;
    //Check if packet is a duplicate
    Packet *p = pktPool->GetPacket(pktID);
    if (p == NULL) {//packet isn't found in the packet pool
        printf("Error: Packet %d doesn't exist in Packet Pool!Aborting..\n", pktID);
        exit(1);
    }
    Header *h = p->getHeader();
    if (h->IsDuplicate() == true) {
        PacketID = h->GetOriginal();
    } else {
        PacketID = pktID;
    }
    //Sanity check
    if (rTime < 0 || p->GetStartTime() < 0) {
        printf("%f:Node %d received new packet with ID:%d and type %d from %d\n", rTime, this->NodeID, p->getID(),
               p->getType(), h->GetprevHop());
        exit(1);
    }
    //Recognize packet type

    switch (p->getType()) {
        case DATA_PACKET: {
            ReceptionData(h, p, pktID, rTime, PacketID);
            break;
        }
        case DIRECT_SUMMARY_PACKET: {
            ReceptionDirectSummary(h, p, pktID, rTime);
            break;
        }
        case DIRECT_REQUEST_PACKET: {
            ReceptionDirectRequest(h, p, pktID, rTime);
            break;
        }
        case REQUEST_CONTACTS_PACKET: {
            ReceptionRequestContacts(h, p, pktID, rTime);
            break;
        }
        case CONTACTS_PACKET: {
            ReceptionContacts(h, p, pktID, rTime);
            break;
        }
        case REQ_RSPM: {
            ReceptionRequestRSPM(h, p, pktID, rTime);
            break;
        }
        case IND_RSPM: {
            ReceptionRSPM(h, p, pktID, rTime);
            break;
        }
        case PKTS_DESTS: {
            ReceptionSummary(h, p, pktID, rTime);
            break;
        }
        case REQUEST_PKTS_UTILS: {
            ReceptionRequest(h, p, pktID, rTime);
            break;
        }
        case REQUEST_PKTS_MULTI_UTILS: {
            if (TwoDOn) {//Added by epap
                ReceptionRequestSimBetTS2D(h, p, pktID, rTime);
            } else {
                ReceptionRequestSimBetTS(h, p, pktID, rTime);
            }
            break;
        }
        case ANTIPACKET: {
            ReceptionAntipacket(h, p, pktID, rTime);
            break;
        }
        case ANTIPACKET_RESPONSE: {
            ReceptionAntipacketResponse(h, p, pktID, rTime);
            break;
        }
        case REQUEST_BUFFER_INFO: {
            ReceptionBufferReq(h, p, pktID, rTime);
            break;
        }
        case BUFFER_INFO: {
            ReceptionBufferRsp(h, p, pktID, rTime);
            break;
        }
        case DEST_PREDICT_PACKET: {
            ReceptionDPT(h, p, pktID, rTime);
            break;
        }
        default: {
            printf("Error: Unknown packet type!CustomCbRDF Forwarding protocol is not using this type of packet..Aborting...\n");
            exit(1);
        }
    }
    return;
}


void CustomCbRDF::SendContactRequest(double CTime, int NID) {
    //Create new request contacts packet
    Packet *ReqPacket = new ReqContacts(CTime, 0);
    Header *Reqh = new BasicHeader(this->NodeID, NID);
    ReqPacket->setHeader(Reqh);
    //Add packet to the packet pool
    pktPool->AddPacket(ReqPacket);
    //Send packet to the new contact
    Mlayer->SendPkt(CTime, this->NodeID, NID, ReqPacket->getSize(), ReqPacket->getID());
#ifdef CustomCbRDF_DEBUG
    printf("%d:Send a request for contacts to node %d\n",this->NodeID,NID);
#endif
    return;
}


void CustomCbRDF::ReceptionRequestContacts(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    if (this->UType == "Prophet") {
        //Create response containing contacts
        double *delivpre = MyDPT->CloneDPT(CurrentTime);
        Packet *Response = new DPs(CurrentTime, 0);
        Response->setContents((void *) delivpre);
        Header *reshd = new BasicHeader(this->NodeID, hd->GetprevHop());
        Response->setHeader(reshd);
        //Add packet to the packet pool
        pktPool->AddPacket(Response);
        //Send packet
        Mlayer->SendPkt(CurrentTime, this->NodeID, hd->GetprevHop(), Response->getSize(), Response->getID());
#ifdef CustomCbRDF_DEBUG
        printf("%f:Node %d send its DPT with ID:%d to %d\n",CurrentTime,this->NodeID,Response->getID(),hd->GetprevHop());
        printf("DPT contents(%d):\n",Set->getNN());
        for(int i=0;i<Set->getNN();i++)
        {
            printf("%d->%f ",i,delivpre[i]);
        }
        printf("\n");
#endif
    } else {
        //request for contacts packet
#ifdef CustomCbRDF_DEBUG
        printf("%f:Node %d received a request for contacts with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
#endif
        //Create response containing contacts
        int *con = Adja->GetMyContacts();
        Packet *Response = new Contacts(CurrentTime, 0);
        Response->setContents((void *) con);
        Header *reshd = new BasicHeader(this->NodeID, hd->GetprevHop());
        Response->setHeader(reshd);
        //Add packet to the packet pool
        pktPool->AddPacket(Response);
        //Send packet
        Mlayer->SendPkt(CurrentTime, this->NodeID, hd->GetprevHop(), Response->getSize(), Response->getID());
#ifdef CustomCbRDF_DEBUG
        printf("%f:Node %d send its contacts with ID:%d to %d\n",CurrentTime,this->NodeID,Response->getID(),hd->GetprevHop());
        printf("Contact contents(%d):\n",con[0]);
        for(int i=1;i<=con[0];i++)
        {
            printf("%d ",con[i]);
        }
        printf("\n");
#endif
    }
    //Delete packet to free memory
    pktPool->ErasePacket(PID);
    return;
}


void CustomCbRDF::ReceptionContacts(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //Packet with contacts
    //Get packet contents
    int *info = (int *) pkt->getContents();
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a packet containing contacts with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
    printf("Contacts contents(%d):\n",info[0]);
    for(int i=1;i<=info[0];i++)
    {
        printf("%d ",info[i]);
    }
    printf("\n");
#endif
    //Update Ego matrix
    for (int i = 1; i <= info[0]; i++) {
        Adja->SetConnection(hd->GetprevHop(), info[i], CurrentTime);
    }
    //Update metrics
    Adja->UpdateAll();
    //Reply with a packet summary
    SendSummary(CurrentTime, hd->GetprevHop());
    //Delete packet to free memory
    pktPool->ErasePacket(PID);
    return;
}


/* ReceptionDPT
 * ------------
 * This method is called when a packet, that contains the delivery probabilities of the node in
 * contact, is received. The proper updates in the receiving node's delivery probability table
 * take place.
 */
void CustomCbRDF::ReceptionDPT(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //packet containing DPT
    //Get packet contents
    double *info = (double *) pkt->getContents();
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a packet containing a DPT with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
    printf("DPT contents(%d):\n",Set->getNN());
    for(int i=0;i<Set->getNN();i++)
    {
        printf("%d->%f ",i,info[i]);
    }
    printf("\n");
#endif
    MyDPT->UpdateDPT(info, hd->GetprevHop(), CurrentTime);
    //Reply with a packet summary
    SendSummary(CurrentTime, hd->GetprevHop());
    //Delete packet to free memory
    pktPool->ErasePacket(PID);
    return;
}


void CustomCbRDF::SendSummary(double CTime, int NID) {
    //Prepare the summary vector for all other packets
//  	int *OtherSummary=Buf->getAllPackets();
    int *OtherSummary = Buf->getPacketsNotDestinedTo(NID);

    struct PktDest *Sum = (struct PktDest *) malloc(sizeof(struct PktDest) * OtherSummary[0]);
    for (int i = 1; i <= OtherSummary[0]; i++) {
        Sum[i - 1].PID = OtherSummary[i];
        Sum[i - 1].Dest = Buf->GetPktDestination(OtherSummary[i]);
    }
    //Create new summary packet
    Packet *SumPacket = new PktDests(CTime, 0);
    SumPacket->setContents((void *) Sum);
    Header *SumHeader = new BasicHeader(this->NodeID, NID);
    SumPacket->setHeader(SumHeader);
    ((PktDests *) SumPacket)->setPktNum(OtherSummary[0]);
    //Add packet to the packet pool
    pktPool->AddPacket(SumPacket);
    //Send packet to the new contact
    Mlayer->SendPkt(CTime, this->NodeID, NID, SumPacket->getSize(), SumPacket->getID());
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d generated new summary packet with ID:%d for Node %d\n",CTime,this->NodeID,SumPacket->getID(),NID);
    printf("Summary contents(%d):\n",OtherSummary[0]);
    for(int i=0;i<OtherSummary[0];i++)
    {
        printf("%d(Dest:%d) ",Sum[i].PID,Sum[i].Dest);
    }
    printf("\n");
#endif
    free(OtherSummary);
    return;
}


void CustomCbRDF::ReceptionSummary(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //Summary packet
    //Get packet contents
    struct PktDest *summary = (struct PktDest *) pkt->getContents();
    int numContents = ((PktDests *) pkt)->getPktNum();
    int *req = (int *) malloc(sizeof(int) * numContents);
    int *reqDest = (int *) malloc(sizeof(int) * numContents);
    bool *haveIt = (bool *) malloc(sizeof(bool) * numContents);
    for (int i = 0; i < numContents; i++) {
        req[i] = summary[i].PID;
        reqDest[i] = summary[i].Dest;
        if (Buf->PacketExists(summary[i].PID)) {
            haveIt[i] = true;
        } else {//Packet does not exist in node's buffer
            haveIt[i] = false;
        }
    }
    if (UType == "SimBetTS" || UType == "SimBet") {
        prepareReqSimBetTS(CurrentTime, hd->GetprevHop(), numContents, req, reqDest, haveIt);
    } else {
        prepareReq(CurrentTime, hd->GetprevHop(), numContents, req, reqDest, haveIt);
    }
    //free memory
    free(req);
    free(reqDest);
    free(haveIt);
    //Delete summary packet
    pktPool->ErasePacket(PID);
}


void
CustomCbRDF::prepareReqSimBetTS(double CurrentTime, int encID, int numContents, int *req, int *reqDest, bool *haveIt) {
    //initialize the data structure that will form the reply
    struct PktMultiUtil *RList = (struct PktMultiUtil *) malloc(sizeof(struct PktMultiUtil) * numContents);
    //for each packet set
    for (int i = 0; i < numContents; i++) {
        //packet ID, whether this node owns it or not
        RList[i].PID = req[i];
        RList[i].Exists = haveIt[i];
        RList[i].mode = 0;
        double Bt = Adja->getBet();
        double Sm = Adja->getSim(reqDest[i]);
        double Freq = Adja->getFreq(reqDest[i]);
        double Intim = Adja->getIntimacy(reqDest[i]);
        double Rec = Adja->getRecency(reqDest[i], CurrentTime);

        RList[i].Sim = (Sm);
        RList[i].Bet = (Bt);
        RList[i].Frequency = (Freq);
        RList[i].Intimacy = (Intim);
        RList[i].Recency = (Rec);

    }
    //Create a packet request as a response
    Packet *ReqPacket = new PktMultiUtils(CurrentTime, 0);
    Header *head = new BasicHeader(this->NodeID, encID);
    ReqPacket->setHeader(head);
    ReqPacket->setContents((void *) RList);
    ((PktMultiUtils *) ReqPacket)->setPktNum(numContents);
    //Add packet to the packet pool
    pktPool->AddPacket(ReqPacket);
    //Send packet to the new contact
    Mlayer->SendPkt(CurrentTime, this->NodeID, encID, ReqPacket->getSize(), ReqPacket->getID());
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a summary packet with ID:%d from %d\n",CurrentTime,this->NodeID,PID,encID);
    printf("Summary contents(%d):\n",numContents);
    for(int i=0;i<numContents;i++)
    {
        printf("%d(Dest:%d) ",summary[i].PID,summary[i].Dest);
    }
    printf("\n");
    fflush(stdout);
#endif
    return;
}

void CustomCbRDF::prepareReq(double CurrentTime, int encID, int numContents, int *req, int *reqDest, bool *haveIt) {
    double MyBet = 0.0;
    int locatedKmobject;
    if (Adja) {
        if (this->UType == "Bet") {
            MyBet = Adja->getBet();
        }
    }

    //---------------------------------------------Added by giorg---------------------------------------------//
    //CustomPktUtil(Packet.h) domh gia na apothikeuei kai tis 3 metrikes poy xreiazomaste na steiloyme

    //struct PktUtil *RList=(struct PktUtil *)malloc(sizeof(struct PktUtil)*numContents);

    struct CustomPktUtil *RList = (struct CustomPktUtil *) malloc(sizeof(struct CustomPktUtil) * numContents);
    double util = 0.0;
    for (int i = 0; i < numContents; i++) {
        RList[i].PID = req[i];
        RList[i].Exists = haveIt[i];
        RList[i].mode = 0;
        if (MyDPT && this->UType == "Prophet") {
            util = MyDPT->getDPto(reqDest[i], CurrentTime);
        } else if (Adja) {
            if (this->UType == "Bet") {
                util = MyBet;
            } else if (this->UType == "Sim") {
                util = Adja->getSim(reqDest[i]);
            } else {
                util = 0.0;
            }
        } else {
            util = Util->get(reqDest[i], CurrentTime);
        }
        RList[i].Util = util;

//---------------------------------------------Added by giorg---------------------------------------------//

//elegxos gia to an yparxei to kmeans object gia to sigkerimeno dest allios arxikopoihsh me arnhtikes times//

        RList[i].avr = -100;
        RList[i].sampleNumber = 0;
        RList[i].standDeviation =-100;

        for (int km = 0; km < (int) KmRouting.size(); km++) {
            if (KmRouting[km].getDest() == reqDest[i]) {
                locatedKmobject = km;
                RList[i].avr = KmRouting[locatedKmobject].average();
                RList[i].sampleNumber = KmRouting[locatedKmobject].getElementcnt();
                RList[i].standDeviation = KmRouting[locatedKmobject].stdev();
                break;
            }

        }
//-----------------------------------------------------------------------------------------------------//




    }
    //Create a packet request as a response
    Packet *ReqPacket = new PktUtils(CurrentTime, 0);
    Header *head = new BasicHeader(this->NodeID, encID);
    ReqPacket->setHeader(head);
    ReqPacket->setContents((void *) RList);
    ((PktUtils *) ReqPacket)->setPktNum(numContents);
    //clean this - not needed
    ((PktUtils *) ReqPacket)->setDIutil(0.0);
    //Add packet to the packet pool
    pktPool->AddPacket(ReqPacket);
    //Send packet to the new contact
    Mlayer->SendPkt(CurrentTime, this->NodeID, encID, ReqPacket->getSize(), ReqPacket->getID());
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a summary packet with ID:%d from %d\n",CurrentTime,this->NodeID,PID,encID);
    printf("Summary contents(%d):\n",numContents);
    for(int i=0;i<numContents;i++)
    {
        printf("%d(Dest:%d) ",summary[i].PID,summary[i].Dest);
    }
    printf("\n");
#endif
    return;
}


void CustomCbRDF::ReceptionRequest(Header *hd, Packet *pkt, int PID, double CurrentTime) {

    //Received packet containing pkts requested and their utilities
    struct CustomPktUtil *RqList = (struct CustomPktUtil *) pkt->getContents();

    int PktNum = ((PktUtils *) pkt)->getPktNum();
    //split packet contents to three arrays
    int *Requests = (int *) malloc(sizeof(int) * (PktNum + 1));
    double *Utils = (double *) malloc(sizeof(double) * (PktNum + 1));
    bool *OwnedByOther = (bool *) malloc(sizeof(bool) * (PktNum + 1));

    //---------------------------------------------Added by giorg---------------------------------------------//
    // prostheto three arrays gia na apothikeuso ta extra data poy periexei to Received packet.
    double *avrPrice = (double *) malloc(sizeof(double)*(PktNum + 1));
    double *stdDev = (double *) malloc(sizeof(double)*(PktNum + 1));
    int *sample = (int *) malloc(sizeof(int) * (PktNum + 1));

    //-----------------------------------------------------------------------------------------------------//

    //---------------------------------------------Added by Vag---------------------------------------------//
    int *pktsdests = (int *) malloc(sizeof(int) * (PktNum + 1));
    int *pktsdecisions = (int *) malloc(sizeof(int) * (PktNum +
                                                       1)); //Store decision for each packet: -1->training, 0->do not copy(not in training), 1->copy (not in training)

    bool destinationAlreadyAdded = false;
    int locatedKmobject = -1;
    bool multiplepacketsperdst = false;
    //double clusterofthreshold=0.0;
    //double clusterofother=0.0;
    int clusterofthresholdrank = 0;
    int clusterofotherrank = 0;
    int clusterofmerank = 0;
    pktsdecisions[0] = PktNum;
    //-----------------------------------------------------------------------------------------------------//

    Requests[0] = PktNum;
    Utils[0] = PktNum;
    OwnedByOther[0] = false;
    //---------------------------------------------Added by giorg---------------------------------------------//
    avrPrice[0]=PktNum;
    stdDev[0]=PktNum;
    sample[0]=PktNum;
    //-----------------------------------------------------------------------------------------------------//

    for (int i = 1; i <= PktNum; i++) {
        Requests[i] = RqList[i - 1].PID;
        Utils[i] = RqList[i - 1].Util;
        OwnedByOther[i] = RqList[i - 1].Exists;

        //---------------------------------------------Added by giorg---------------------------------------------//
        avrPrice[i] = RqList[i - 1].avr;
        sample[i]=RqList[i - 1].sampleNumber;
        stdDev[i] =RqList[i - 1].standDeviation;
        //-----------------------------------------------------------------------------------------------------//


        //---------------------------------------------Added by Vag---------------------------------------------//
        pktsdecisions[i] = -1;
        pktsdests[i] = -1;

        if (Utils[i] != Utils[i])
            continue; //This is a hack for the case that Utils[i] is not a number (appears when using Cambridge trace). In this case the corresponding packets is not forwarded.

        clusterofthresholdrank = 0;
        clusterofotherrank = 0;
        clusterofmerank = 0;
        pktsdests[i] = Buf->getPacketDest(Requests[i]);

        locatedKmobject = -1;

        if (DstDepend) {
            multiplepacketsperdst = false;
            for (int jjj = 1; jjj <= i - 1; jjj++) {
                if (pktsdests[i] == pktsdests[jjj]) {
                    multiplepacketsperdst = true;
                    break;
                }
            }
            destinationAlreadyAdded = false;
            //Locate Km object for destination or create it if it does not exist

            for (int km = 0; km < (int) KmRouting.size(); km++) {
                if (KmRouting[km].getDest() == pktsdests[i]) {
                    destinationAlreadyAdded = true;
                    locatedKmobject = km;
                    break;
                }
            }
            if (destinationAlreadyAdded == false) {
                locatedKmobject = (int) KmRouting.size();
                KmRouting.push_back(Kmeans(this->NodeID, pktsdests[i], isKmeansnormalized, isLVQon, KmeansUPeriod,
                                           isKmeansweighted));
            }
        } else {
            if ((int) KmRouting.size() == 0)
                KmRouting.push_back(
                        Kmeans(this->NodeID, -1, isKmeansnormalized, isLVQon, KmeansUPeriod, isKmeansweighted));
            locatedKmobject = 0;
        }

        //Either store new utility for destination or update based on LVQ
          //-----------------------------------------------------------Added by giorg------------------------------------------------------//


      /*  if(KmRouting[locatedKmobject].isInTraining()){
            if (((!multiplepacketsperdst)&&(DstDepend)) || ((DstDepend == false)&&(i == 1))) KmRouting[locatedKmobject].setElementOfUtilities(Utils[i]);
        } else {
            if (((!multiplepacketsperdst)&&(DstDepend)) || ((DstDepend == false)&&(i == 1))) KmRouting[locatedKmobject].UpdatewithNewElement(Utils[i]);
        }*/

        //metatroph ton SetElementOfUtilities, UpdateWithNewElement etsi oste na dexontai ta 3 data(avrPrice , stdDev , sample)

        if (KmRouting[locatedKmobject].isInTraining()) {
            if (((!multiplepacketsperdst) && (DstDepend)) || ((DstDepend == false) && (i == 1)))
                KmRouting[locatedKmobject].CustomSetElementOfUtilities(Utils[i],avrPrice[i],stdDev[i],sample[i]);
        } else {
            if (((!multiplepacketsperdst) && (DstDepend)) || ((DstDepend == false) && (i == 1)))
                KmRouting[locatedKmobject].CustomUpdatewithNewElement(Utils[i],avrPrice[i],stdDev[i],sample[i]);
        }

       //----------------------------------------------------------------------------------------------------------------------------------//
        //printf("Node with id %d for destination: %d have size: %d, avgUtil: %f, stdDev: %f.\n",this->NodeID,pktsdests[i],KmRouting[locatedKmobject].getElementcnt(),KmRouting[locatedKmobject].average(),KmRouting[locatedKmobject].stdev());

        if (KmRouting[locatedKmobject].isInTraining() == false) {
            clusterofthresholdrank = KmRouting[locatedKmobject].returnClusterRank(
                    (Buf->getPacketData(Requests[i]))->GetMaxUtil());
            clusterofotherrank = KmRouting[locatedKmobject].returnClusterRank(Utils[i]);
            clusterofmerank = KmRouting[locatedKmobject].returnClusterRank(getmyUtil(pktsdests[i], CurrentTime));

            //1 Best strategy so far: forward only to a better quality group but boost (by replicating to peers) in the first hop.
            //if ((clusterofotherrank < clusterofthresholdrank) || ((clusterofotherrank == clusterofthresholdrank)&&(Buf->GetHops(Requests[i])==0))) {
            //2 Second preferable strategy: forward only to a better quality group but boost (by replicating to peers) in the first k hops. Note: Try to define k based on the total number of clusters and the current cluster.
            //if ((clusterofotherrank < clusterofthresholdrank) || ((clusterofotherrank == clusterofthresholdrank)&&(Buf->GetHops(Requests[i])<2))) {
            //3 Probably better solution: Very promising. Smaller overhead reduction but smaller delay compared to other solutions
            //if ((clusterofotherrank < clusterofthresholdrank) || ((clusterofotherrank == clusterofthresholdrank)&&((clusterofmerank == clusterofthresholdrank)||(Buf->GetHops(Requests[i]) == 0)))) {
            //**************************************************
            //NOTE: my current rank may be better than the threshold rank
            //if (clusterofmerank < clusterofthresholdrank) printf("Very strange.\n");
            //**************************************************
            //4 Another promising one
            //if ((clusterofotherrank < clusterofthresholdrank) || ((clusterofotherrank == clusterofthresholdrank)&&((clusterofmerank == clusterofthresholdrank)||(clusterofmerank == 1)))) {
            //5 Testing now
            //Probably good because it is a mixture of policies (3) and (4)
            if ((clusterofotherrank < clusterofthresholdrank) || ((clusterofotherrank == clusterofthresholdrank) &&
                                                                  ((clusterofmerank == clusterofthresholdrank) ||
                                                                   (clusterofmerank ==
                                                                    1)/*||(Buf->GetHops(Requests[i]) == 0)*/))) {
                pktsdecisions[i] = 1;
            } else {
                pktsdecisions[i] = 0;
            }
        }
        //-----------------------------------------------------------------------------------------------------//
    }
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a request packet with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
    printf("Request contents(%d):\n",PktNum);
    for(int i=0;i<PktNum;i++)
    {
        printf("%d(Util:%f,",Requests[i+1],Utils[i+1]);
        if(OwnedByOther[i+1])
        {
            printf("E) ");
        }
        else
        {
            printf("NoE) ");
        }
    }
    printf("\n");
    printf("Buffer contents\n");
    printf("---------------\n");
    Buf->PrintPkts();
#endif
    //---------------------------------------------Added by Vag---------------------------------------------//
    int *pktsToSend = Buf->getPacketsLowUtilandKMeans(Requests, Utils, OwnedByOther, pktsdecisions);
    //-----------------------------------------------------------------------------------------------------//
#ifdef CustomCbRDF_DEBUG
    printf("Update completed!\n");
    printf("Buffer contents\n");
    printf("---------------\n");
    Buf->PrintPkts();
#endif
    //Schedule packets for transmission
    for (int i = 1; i <= pktsToSend[0]; i++) {
        sch->addPacket(pktsToSend[i], NULL);
    }
    free(pktsToSend);
    //Apply scheduling
    int *outgoing = sch->getOutgoingPackets();
    //Apply congestion control
    outgoing = CC->filterPackets(outgoing);
    if (outgoing) {
        for (int i = 1; i <= outgoing[0]; i++) {
            SendPacket(CurrentTime, outgoing[i], hd->GetprevHop(), 1);
        }
        free(outgoing);
    }
    //free memory
    free(Requests);
    free(Utils);
    free(OwnedByOther);

    //---------------------------------------------Added by Vag---------------------------------------------//
    free(pktsdests);
    free(pktsdecisions);
    //-----------------------------------------------------------------------------------------------------//


//---------------------------------------------Added by giorg---------------------------------------------//
    free(avrPrice);
    free(stdDev);
    free(sample);
//-----------------------------------------------------------------------------------------------------//


    pktPool->ErasePacket(PID);
    return;
}

//-----------------------------------------------Added by epap-----------------------------------//
/* Implementation using SimBet and SimBetTS as defined by the respective protocols - Problem with normalization*/
void CustomCbRDF::ReceptionRequestSimBetTS(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //Received packet containing pkts requested and their utilities
    struct PktMultiUtil *RqList = (struct PktMultiUtil *) pkt->getContents();
    int PktNum = ((PktMultiUtils *) pkt)->getPktNum();

    //------------------------Added by epap-------------------------------------------//
    int *pktsdests = (int *) malloc(sizeof(int) * (PktNum + 1));
    int *pktsdecisions = (int *) malloc(sizeof(int) * (PktNum +
                                                       1)); //Store decision for each packet: -1->training, 0->do not copy(not in training), 1->copy (not in training)

    bool destinationAlreadyAdded = false;
    int locatedKmobject = -1;
    bool multiplepacketsperdst = false;

    int clusterofthresholdrankA = 0;
    int clusterofthresholdrankB = 0;
    int clusterofotherrank = 0;
    int clusterofmerank = 0;
    pktsdecisions[0] = PktNum;

    double readUtilreltothresh = -1000.0;
    double threshUtilreltocontact = -1000.0;
    double threshUtilreltome = -1000.0;
    double myUtilreltothresh = -1000.0;

    for (int i = 1; i <= PktNum; i++) {

        pktsdecisions[i] = -1;
        pktsdests[i] = -1;

        pktsdests[i] = Buf->getPacketDest(RqList[i - 1].PID);
        if (this->UType == "SimBet") {
            threshUtilreltocontact = Buf->CalculateSimBetUtility(
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Sim"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Bet"), RqList[i - 1].Sim, RqList[i - 1].Bet);
            readUtilreltothresh = Buf->CalculateSimBetUtility(RqList[i - 1].Sim, RqList[i - 1].Bet,
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Sim"),
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Bet"));
            myUtilreltothresh = Buf->CalculateSimBetUtility(Adja->getSim(pktsdests[i]), Adja->getBet(),
                                                            (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Sim"),
                                                            (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Bet"));
            threshUtilreltome = Buf->CalculateSimBetUtility((Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Sim"),
                                                            (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Bet"),
                                                            Adja->getSim(pktsdests[i]), Adja->getBet());
        } else {
            threshUtilreltocontact = Buf->CalculateSimBetTSUtility(
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Sim"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Bet"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Freq"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Intimacy"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Recency"), RqList[i - 1].Sim,
                    RqList[i - 1].Bet, RqList[i - 1].Frequency, RqList[i - 1].Intimacy, RqList[i - 1].Recency);
            readUtilreltothresh = Buf->CalculateSimBetTSUtility(RqList[i - 1].Sim, RqList[i - 1].Bet,
                                                                RqList[i - 1].Frequency, RqList[i - 1].Intimacy,
                                                                RqList[i - 1].Recency,
                                                                (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                        "Sim"),
                                                                (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                        "Bet"),
                                                                (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                        "Freq"),
                                                                (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                        "Intimacy"),
                                                                (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                        "Recency"));
            myUtilreltothresh = Buf->CalculateSimBetTSUtility(Adja->getSim(pktsdests[i]), Adja->getBet(),
                                                              Adja->getFreq(pktsdests[i]),
                                                              Adja->getIntimacy(pktsdests[i]),
                                                              Adja->getRecency(pktsdests[i], CurrentTime),
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Sim"),
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Bet"),
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Freq"),
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Intimacy"),
                                                              (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil(
                                                                      "Recency"));
            threshUtilreltome = Buf->CalculateSimBetTSUtility(
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Sim"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Bet"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Freq"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Intimacy"),
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Recency"), Adja->getSim(pktsdests[i]),
                    Adja->getBet(), Adja->getFreq(pktsdests[i]), Adja->getIntimacy(pktsdests[i]),
                    Adja->getRecency(pktsdests[i], CurrentTime));
        }
        if (readUtilreltothresh != readUtilreltothresh)
            continue; //This is a hack for the case that Utils[i] is not a number (appears when using Cambridge trace). In this case the corresponding packets is not forwarded.


        clusterofthresholdrankA = 0;
        clusterofthresholdrankB = 0;
        clusterofotherrank = 0;
        clusterofmerank = 0;

        locatedKmobject = -1;

        multiplepacketsperdst = false;
        for (int jjj = 1; jjj <= i - 1; jjj++) {
            if (pktsdests[i] == pktsdests[jjj]) {
                multiplepacketsperdst = true;
                break;
            }
        }
        destinationAlreadyAdded = false;
        //Locate Km object for destination or create it if it does not exist
        for (int km = 0; km < (int) KmRouting.size(); km++) {
            if (KmRouting[km].getDest() == pktsdests[i]) {
                destinationAlreadyAdded = true;
                locatedKmobject = km;
                break;
            }
        }
        if (destinationAlreadyAdded == false) {
            locatedKmobject = (int) KmRouting.size();
            KmRouting.push_back(
                    Kmeans(this->NodeID, pktsdests[i], isKmeansnormalized, isLVQon, KmeansUPeriod, isKmeansweighted));
        }
        //Either store new utility for destination or update
        if (KmRouting[locatedKmobject].isInTraining()) {
            if (!multiplepacketsperdst) KmRouting[locatedKmobject].setElementOfUtilities(readUtilreltothresh);
        } else {
            if (!multiplepacketsperdst) KmRouting[locatedKmobject].UpdatewithNewElement(readUtilreltothresh);
        }

        if (KmRouting[locatedKmobject].isInTraining() == false) {
            clusterofthresholdrankA = KmRouting[locatedKmobject].returnClusterRank(threshUtilreltocontact);
            clusterofthresholdrankB = KmRouting[locatedKmobject].returnClusterRank(threshUtilreltome);
            clusterofotherrank = KmRouting[locatedKmobject].returnClusterRank(readUtilreltothresh);
            clusterofmerank = KmRouting[locatedKmobject].returnClusterRank(myUtilreltothresh);

            if ((clusterofotherrank < clusterofthresholdrankA) || ((clusterofotherrank == clusterofthresholdrankA) &&
                                                                   ((clusterofmerank == clusterofthresholdrankB) ||
                                                                    (clusterofmerank == 1)))) {
                pktsdecisions[i] = 1;
            } else {
                pktsdecisions[i] = 0;
            }
        }

    }

    //--------------------------------------------------------------------------------//
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a request packet with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
    printf("Request contents(%d):\n",PktNum);
    for(int i=0;i<PktNum;i++)
    {
        printf("%d(Sim:%f,Bet:%f,Freq:%f,Intimacy:%f,Recency:%f) ",RqList[i].PID,RqList[i].Sim,RqList[i].Bet,RqList[i].Frequency,RqList[i].Intimacy,RqList[i].Recency);
    }
    printf("\n");
    fflush(stdout);
#endif
    int *pktsToSend = NULL;
    //----------------------------------Modified by epap------------------------------//
    if (this->UType == "SimBet") {
        pktsToSend = Buf->getPacketsLowUtilSimBetandKMeans(PktNum, RqList, pktsdecisions);
    } else {
        pktsToSend = Buf->getPacketsLowUtilSimBetTSandKMeans(PktNum, RqList, pktsdecisions);
    }
    //--------------------------------------------------------------------------------//
#ifdef CustomCbRDF_DEBUG
    printf("Update completed!\n");
    printf("Buffer contents\n");
    printf("---------------\n");
    Buf->PrintPkts();
#endif
    //Schedule packets for transmission
    for (int i = 1; i <= pktsToSend[0]; i++) {
        sch->addPacket(pktsToSend[i], NULL);
    }
    free(pktsToSend);
    //Apply scheduling
    int *outgoing = sch->getOutgoingPackets();
    //Apply congestion control
    outgoing = CC->filterPackets(outgoing);
    if (outgoing) {
        for (int i = 1; i <= outgoing[0]; i++) {
            SendPacket(CurrentTime, outgoing[i], hd->GetprevHop(), 1);
        }
        free(outgoing);
    }
    //free memory
    //---------------Added by epap---------------//
    free(pktsdests);
    free(pktsdecisions);
    //-------------------------------------------//
    pktPool->ErasePacket(PID);
    return;
}

//Implementation for Kmeans with two metrics
void CustomCbRDF::ReceptionRequestSimBetTS2D(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //Received packet containing pkts requested and their utilities
    struct PktMultiUtil *RqList = (struct PktMultiUtil *) pkt->getContents();
    int PktNum = ((PktMultiUtils *) pkt)->getPktNum();
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a request packet with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
    printf("Request contents(%d):\n",PktNum);
    for(int i=0;i<PktNum;i++)
    {
        printf("%d(Sim:%f,Bet:%f,Freq:%f,Intimacy:%f,Recency:%f) ",RqList[i].PID,RqList[i].Sim,RqList[i].Bet,RqList[i].Frequency,RqList[i].Intimacy,RqList[i].Recency);
    }
    printf("\n");
    fflush(stdout);
#endif

    //---------------------------------------------Added by epap--------------------------------------------//
    int *pktsdests = (int *) malloc(sizeof(int) * (PktNum + 1));
    int *pktsdecisions = (int *) malloc(sizeof(int) * (PktNum +
                                                       1)); //Store decision for each packet: -1->training, 1->contact in better Sim group, 2->contact in better Bet group, 0->contact not in better Sim or Bet group
    pktsdecisions[0] = PktNum;

    int locatedKmobject = -1;
    bool destinationAlreadyAdded = false;
    bool multiplepacketsperdst = false;
    int clusterofthresholdrank = 0;
    int clusterofotherrank = 0;
    int clusterofmerank = 0;
    int bclusterofthresholdrank = 0;
    int bclusterofotherrank = 0;
    int bclusterofmerank = 0;

    if (((int) KmRouting.size() == 0) && (PktNum != 0)) {
        KmRouting.push_back(Kmeans(this->NodeID, -1, isKmeansnormalized, isLVQon, KmeansUPeriod, isKmeansweighted));
    }

    for (int i = 1; i <= PktNum; i++) {
        pktsdecisions[i] = -1;
        pktsdests[i] = -1;

        pktsdests[i] = Buf->getPacketDest(RqList[i - 1].PID);

        clusterofthresholdrank = 0;
        clusterofotherrank = 0;
        clusterofmerank = 0;
        bclusterofthresholdrank = 0;
        bclusterofotherrank = 0;
        bclusterofmerank = 0;

        multiplepacketsperdst = false;
        for (int jjj = 1; jjj <= i - 1; jjj++) {
            if (pktsdests[i] == pktsdests[jjj]) {
                multiplepacketsperdst = true;
                break;
            }
        }
        destinationAlreadyAdded = false;
        locatedKmobject = -1;

        //Locate Km object for destination or create it if it does not exist
        for (int km = 0; km < (int) KmRouting.size(); km++) {
            if (KmRouting[km].getDest() == pktsdests[i]) {
                destinationAlreadyAdded = true;
                locatedKmobject = km;
                break;
            }
        }
        if (destinationAlreadyAdded == false) {
            //if (multiplepacketsperdst) printf("BIG PROBLEM.\n");
            locatedKmobject = (int) KmRouting.size();
            KmRouting.push_back(
                    Kmeans(this->NodeID, pktsdests[i], isKmeansnormalized, isLVQon, KmeansUPeriod, isKmeansweighted));
        }

        //Either store new utility for destination or update based on LVQ
        if (KmRouting[locatedKmobject].isInTraining()) {
            if (!multiplepacketsperdst) KmRouting[locatedKmobject].setElementOfUtilities(RqList[i - 1].Sim);
        } else {
            if (!multiplepacketsperdst) KmRouting[locatedKmobject].UpdatewithNewElement(RqList[i - 1].Sim);
        }

        if ((i == 1) && ((int) KmRouting.size() > 0)) {
            if (KmRouting[0].isInTraining()) {
                KmRouting[0].setElementOfUtilities(RqList[i - 1].Bet);
            } else {
                KmRouting[0].UpdatewithNewElement(RqList[i - 1].Bet);
            }
        }

        if ((KmRouting[locatedKmobject].isInTraining() == false) && (KmRouting[0].isInTraining() == false)) {

            clusterofthresholdrank = (double) KmRouting[locatedKmobject].returnClusterRank(
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Sim"));
            clusterofotherrank = (double) KmRouting[locatedKmobject].returnClusterRank(RqList[i - 1].Sim);
            clusterofmerank = (double) KmRouting[locatedKmobject].returnClusterRank(
                    getmyUtil(pktsdests[i], CurrentTime, 0));

            bclusterofthresholdrank = (double) KmRouting[0].returnClusterRank(
                    (Buf->getPacketData(RqList[i - 1].PID))->GetMaxUtil("Bet"));
            bclusterofotherrank = (double) KmRouting[0].returnClusterRank(RqList[i - 1].Bet);
            bclusterofmerank = (double) KmRouting[0].returnClusterRank(getmyUtil(pktsdests[i], CurrentTime, 1));

            pktsdecisions[i] = 0;

            //if (clusterofmerank != 1){
            //if (clusterofthresholdrank !=  1){
            if ((clusterofthresholdrank == clusterofmerank) && (clusterofthresholdrank != 1)) {
                if (bclusterofotherrank < bclusterofthresholdrank) {
                    pktsdecisions[i] = 2;
                } else if ((bclusterofotherrank == bclusterofthresholdrank) &&
                           (bclusterofthresholdrank == bclusterofmerank)) {
                    pktsdecisions[i] = 2;
                }
            } else if ((clusterofmerank == 1) && (clusterofotherrank == 1)) {
                //pktsdecisions[i]=2;
                if (bclusterofotherrank < bclusterofthresholdrank) {
                    pktsdecisions[i] = 2;
                } else if ((bclusterofotherrank == bclusterofthresholdrank) &&
                           (bclusterofthresholdrank == bclusterofmerank)) {
                    pktsdecisions[i] = 2;
                }

            }


        }

    }
    //-----------------------------------------------------------------------------------------------------//

    int *pktsToSend = NULL;
    if (this->UType == "SimBet") {
        //-----------------------------------Added by epap-----------------------------------------------//
        pktsToSend = Buf->getPacketsLowUtilSimBetand2DKMeans(PktNum, RqList, pktsdecisions);
        //-----------------------------------------------------------------------------------------------//
    } else {
        pktsToSend = Buf->getPacketsLowUtilSimBetTS(PktNum, RqList);
    }
#ifdef CustomCbRDF_DEBUG
    printf("Update completed!\n");
    printf("Buffer contents\n");
    printf("---------------\n");
    Buf->PrintPkts();
#endif
    //Schedule packets for transmission
    for (int i = 1; i <= pktsToSend[0]; i++) {
        sch->addPacket(pktsToSend[i], NULL);
    }
    free(pktsToSend);
    //Apply scheduling
    int *outgoing = sch->getOutgoingPackets();
    //Apply congestion control
    outgoing = CC->filterPackets(outgoing);
    if (outgoing) {
        for (int i = 1; i <= outgoing[0]; i++) {
            SendPacket(CurrentTime, outgoing[i], hd->GetprevHop(), 1);
        }
        free(outgoing);
    }
    //free memory
    pktPool->ErasePacket(PID);
    //---------------------------------------------Added by Vag---------------------------------------------//
    free(pktsdests);
    free(pktsdecisions);
    //-----------------------------------------------------------------------------------------------------//
    return;
}

//-----------------------------------------------------------------------------------------------//


void CustomCbRDF::ReceptionData(Header *hd, Packet *pkt, int PID, double CurrentTime, int RealID) {
    //Data packet
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received new data packet with ID:%d from %d\n",CurrentTime,this->NodeID,RealID,hd->GetprevHop());
#endif
    //case I am the packet originator (packet comes from application layer)
    if (hd->GetSource() == this->NodeID && hd->GetNextHop() == -1) {
        if (this->UType == "SimBetTS" || this->UType == "SimBet") {
            struct SimBetTSmetrics *U = (struct SimBetTSmetrics *) malloc(sizeof(struct SimBetTSmetrics));
            U->Similarity = Adja->getSim(hd->GetDestination());
            U->Betweenness = Adja->getBet();
            U->Frequency = Adja->getFreq(hd->GetDestination());
            U->Intimacy = Adja->getIntimacy(hd->GetDestination());
            U->Recency = Adja->getRecency(hd->GetDestination(), CurrentTime);
            Buf->addPkt(RealID, hd->GetDestination(), hd->GetSource(), CurrentTime, hd->GetHops(), hd->GetprevHop(),
                        pkt->GetStartTime(), U);
        } else {
            Buf->addPkt(RealID, hd->GetDestination(), hd->GetSource(), CurrentTime, hd->GetHops(), hd->GetprevHop(),
                        pkt->GetStartTime());
            double util = 0.0;
            if (MyDPT && this->UType == "Prophet") {
                util = MyDPT->getDPto(hd->GetDestination(), CurrentTime);
            } else if (Adja) {
                if (this->UType == "Bet") {
                    util = Adja->getBet();
                } else if (this->UType == "Sim") {
                    util = Adja->getSim(hd->GetDestination());
                }
                {
                    util = 0.0;
                }
            } else {
                util = Util->get(hd->GetDestination(), CurrentTime);
            }
            (Buf->getPacketData(RealID))->SetMaxUtil(util);
        }
        Stat->pktGen(RealID, hd->GetSource(), hd->GetDestination(), CurrentTime);
        return;
    }
    //If I'm not the next hop
    if (hd->GetNextHop() != this->NodeID) {
        //Update Statistics
        //Stat->incDuplicates();
        if (pkt->AccessPkt() == false) {
            pktPool->ErasePacket(PID);
        }
        return;
    }
    //Check if the destination of the packet is me
    if (hd->GetDestination() == this->NodeID) {
        //The above lines enable checking for duplicates
        if (DM->NoDuplicatesSupport() && DM->isDelivered(RealID)) {
            printf("Problem: Packet %d has been already delivered!\n", RealID);
            exit(1);
        }
        DM->setAsDelivered(RealID);
        //Update Statistics
        Stat->pktRec(hd->GetHops(), (CurrentTime - pkt->GetStartTime()), pkt, pkt->GetStartTime(), false);
        //Garbage Collecting
        if (pkt->AccessPkt() == false) {
            pktPool->ErasePacket(PID);
        }
        return;
    }
    //I am the next hop
    if (Buf->PacketExists(RealID)) {
        printf("[Error]: Node %d received a packet with ID %d from node %d that already exists in its buffer\n",
               this->NodeID, RealID, hd->GetprevHop());
        exit(EXIT_FAILURE);
    } else {
        if (this->UType == "SimBetTS" || this->UType == "SimBet") {
            double localSim = Adja->getSim(hd->GetDestination());
            double localBet = Adja->getBet();
            double localFreq = Adja->getFreq(hd->GetDestination());
            double localIntim = Adja->getIntimacy(hd->GetDestination());
            double localRec = Adja->getRecency(hd->GetDestination(), CurrentTime);
            struct SimBetTSmetrics *U = (struct SimBetTSmetrics *) malloc(sizeof(struct SimBetTSmetrics));
            U->Similarity = localSim;
            U->Betweenness = localBet;
            U->Frequency = localFreq;
            U->Intimacy = localIntim;
            U->Recency = localRec;
            //add packet to the buffer
            Buf->addPkt(RealID, hd->GetDestination(), hd->GetSource(), CurrentTime, hd->GetHops(), hd->GetprevHop(),
                        pkt->GetStartTime(), U);
            Stat->incTimesAsRelayNode(pkt->GetStartTime());
        } else {
            Buf->addPkt(RealID, hd->GetDestination(), hd->GetSource(), CurrentTime, hd->GetHops(), hd->GetprevHop(),
                        pkt->GetStartTime());
            Stat->incTimesAsRelayNode(pkt->GetStartTime());
            double util = 0.0;
            if (MyDPT && this->UType == "Prophet") {
                util = MyDPT->getDPto(hd->GetDestination(), CurrentTime);
            } else if (Adja) {
                if (this->UType == "Bet") {
                    util = Adja->getBet();
                } else if (this->UType == "Sim") {
                    util = Adja->getSim(hd->GetDestination());
                } else {
                    util = 0.0;
                }
            } else {
                util = Util->get(hd->GetDestination(), CurrentTime);
            }
            (Buf->getPacketData(RealID))->SetMaxUtil(util);
        }
    }
    if (pkt->AccessPkt() == false) {
        pktPool->ErasePacket(PID);
    }
    return;
}


void CustomCbRDF::SendPacket(double STime, int pktID, int nHop, int RepValue) {
// 	printf("%d:Sending packet %d to node %d\n",this->NodeID,pktID,nHop);
    Packet *p = pktPool->GetPacket(pktID);
    if (p == NULL) {//packet isn't found in the packet pool
        printf("Error: Packet %d doesn't exist in Packet Pool!Aborting..\n", pktID);
        exit(1);
    }
    //Duplicate the packet first
    Packet *newPkt = p->Duplicate(Buf->GetHops(pktID));
    newPkt->getHeader()->SetNextHop(nHop);
    newPkt->getHeader()->SetprevHop(this->NodeID);
    newPkt->getHeader()->SetRep(RepValue);
    pktPool->AddPacket(newPkt);
    //Then, inform current neighbors about the new packet
    int CurrentN = Mlayer->BroadcastPkt(STime, this->NodeID, newPkt->getSize(), newPkt->getID());
    //Garbage Collecting
    if (CurrentN > 0) {//Set access attribute to safely delete packet later
        newPkt->SetRecipients(CurrentN);
        //Update statistics
        if (newPkt->getHeader()->GetDestination() == nHop) {
            Stat->incHandovers(Buf->GetPktCreationTime(pktID));
        }
        Stat->incForwards(pktID, Buf->GetPktCreationTime(pktID));
    } else {//Cancel Broadcast and delete packet - case there are no neighbors
        pktPool->ErasePacket(newPkt->getID());
    }
    return;
}

void CustomCbRDF::SendRSPMRequest(double CTime, int NID) {
    //Create new request contacts packet
    Packet *ReqPacket = new ReqRSPM(CTime, 0);
    Header *Reqh = new BasicHeader(this->NodeID, NID);
    ReqPacket->setHeader(Reqh);
    //Add packet to the packet pool
    pktPool->AddPacket(ReqPacket);
    //Send packet to the new contact
    Mlayer->SendPkt(CTime, this->NodeID, NID, ReqPacket->getSize(), ReqPacket->getID());
#ifdef CustomCbRDF_DEBUG
    printf("%d:Send a request for RSPM values to node %d\n",this->NodeID,NID);
#endif
    return;
}

void CustomCbRDF::ReceptionRequestRSPM(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //request for RSPM values
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a request for RSPM values with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
#endif
    //Delete packet to free memory
    pktPool->ErasePacket(PID);
    //Create response containing contacts
    double *val = ((SPM *) Util)->getRSPMfor(hd->GetprevHop());
    Packet *Response = new IndRSPM(CurrentTime, 0);
    Response->setContents((void *) val);
    Header *reshd = new BasicHeader(this->NodeID, hd->GetprevHop());
    Response->setHeader(reshd);
    //Add packet to the packet pool
    pktPool->AddPacket(Response);
    //Send packet
    Mlayer->SendPkt(CurrentTime, this->NodeID, hd->GetprevHop(), Response->getSize(), Response->getID());
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d send its RSPM values (regarding node %d) with ID:%d to %d\n",CurrentTime,this->NodeID,hd->GetprevHop(),Response->getID(),hd->GetprevHop());
    printf("RSPM values: ");
    if(val)
    {
        for(int i=0;i<Set->getNN();i++)
        {
            printf("%d(%f) ",i,val[i]);
        }
        printf("\n");
    }
    else
    {
        printf("None\n");
    }
#endif
    return;
}


void CustomCbRDF::ReceptionRSPM(Header *hd, Packet *pkt, int PID, double CurrentTime) {
    //Packet with RSPM values
    //Get packet contents
    double *info = (double *) pkt->getContents();
#ifdef CustomCbRDF_DEBUG
    printf("%f:Node %d received a packet containing RSPM values with ID:%d from %d\n",CurrentTime,this->NodeID,PID,hd->GetprevHop());
    printf("RSPM values: ");
    if(info)
    {
        for(int i=0;i<Set->getNN();i++)
        {
            printf("%d(%f) ",i,info[i]);
        }
        printf("\n");
    }
    else
    {
        printf("None\n");
    }
#endif
    //Update my RSPM values
    ((SPM *) Util)->updateRSPM(hd->GetprevHop(), info);
    //Delete packet to free memory
    pktPool->ErasePacket(PID);
    //Reply with a packet summary
    SendSummary(CurrentTime, hd->GetprevHop());
    return;
}

//------------------------------------------------------Added by Vag-----------------------------------------------------------------//
double CustomCbRDF::getmyUtil(int dest, double CurrentTime) {

    double myUtil;

    if (MyDPT && this->UType == "Prophet") {
        myUtil = MyDPT->getDPto(dest, CurrentTime);
    } else if (Adja) {
        if (this->UType == "Bet") {
            myUtil = Adja->getBet();
        } else if (this->UType == "Sim") {
            myUtil = Adja->getSim(dest);
        } else {
            myUtil = 0.0;
        }
    } else {
        myUtil = Util->get(dest, CurrentTime);
    }

    return myUtil;

}

double CustomCbRDF::getmyUtil(int dest, double CurrentTime, int SSIM) {

    double myUtil = 0.0;

    if (MyDPT && this->UType == "Prophet") {
        myUtil = MyDPT->getDPto(dest, CurrentTime);
    } else if (Adja) {
        if (SSIM == 0) {
            myUtil = Adja->getSim(dest);
        }
        if (SSIM == 1) {
            myUtil = Adja->getBet();
        }
    } else {
        myUtil = Util->get(dest, CurrentTime);
    }

    return myUtil;

}
//-----------------------------------------------------------------------------------------------------------------------------------//