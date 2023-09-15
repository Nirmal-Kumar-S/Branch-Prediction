#ifndef SIM_BP_H
#define SIM_BP_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <inttypes.h>
#include <math.h>
#include <set>
#include <iomanip>
#include <cstring>
#include "sim_bp.h"

using namespace std;

typedef struct bp_params
{
    unsigned long int K;
    unsigned long int M1;
    unsigned long int M2;
    unsigned long int N;
    char* bp_name;
}bp_params;

class predictor
{
    private:
    int gbh;
    int b_index, g_index, h_index;
    int M1,M2,N,K;
    int mode;
    int *bimodal = NULL;
    int *gshare  = NULL;
    int *hybrid = NULL;

    public:
    int mispredictions;

    predictor(int P_M1, int P_M2, int P_N, int P_K) // Initialize all values of a branch predictor
    {
        gbh = 0;
        M1 = P_M1;
        M2 = P_M2;
        N = P_N;
        K = P_K;
        b_index = g_index = h_index = 0;
        mispredictions = 0;
        int size = 0;
        if (P_M2 > 0 && P_K == 0)
        {   // Bimodal
            mode = 0;
            size = pow(2,M2);
            bimodal = new int[size];
            for (int i = 0; i < size; i++)
            {
                bimodal[i] = 2;
            }   
        } 
        else if (P_M1 > 0 && P_N > 0 && P_K == 0)
        {   // Gshare
            mode = 1;
            size = pow(2,M1);
            gshare = new int[size];
            for (int i = 0; i < size; i++)
            {
                gshare[i] = 2;
            }
        }
        else if (P_M1 > 0 && P_M2 > 0 && P_N > 0 && P_K > 0)         
        {   // Hybrid
            mode = 2;
            size = pow(2,M1);
            gshare = new int[size];
            for (int i =0 ; i < size; i++)
            {    
                gshare[i] = 2;
            }
            size = pow(2,M2);
            bimodal = new int[size];
            for (int i = 0; i < size; i++)
            {    
                bimodal[i] = 2;
            }
            size =  pow(2,K);
            hybrid = new int[size];
            for (int i = 0; i < size; i++)
            {
                hybrid[i] = 1;
            }
        }
    }

    void access(unsigned long int addr, int actual_outcome)
    {
        addr = addr>>2;
        int b_outcome = 0, g_outcome = 0;
        calculate_index(addr);
        if(mode == 0)         
        {   //Bimodal
            if (bimodal[b_index] >= 2)
                b_outcome = 1;

            if(actual_outcome == 1 && bimodal[b_index] < 3)
                bimodal[b_index]++;
            else if (actual_outcome == 0 && bimodal[b_index] > 0)
                bimodal[b_index]--;

            if (b_outcome != actual_outcome)
                mispredictions++;
        }
        else if(mode == 1)         
        {   // Gshare
            if (gshare[g_index] >= 2)
                g_outcome = 1;

            if (actual_outcome == 1 && gshare[g_index] < 3)
                gshare[g_index]++;
            else if (actual_outcome == 0 && gshare[g_index] > 0)
                gshare[g_index]--;
            
            gbh = gbh>>1;
            int temp = actual_outcome<<(N-1);
            gbh = gbh | temp;

            if (g_outcome != actual_outcome)
                mispredictions++;
        }
        else if(mode == 2)         
        {   //Hybrid
            if (bimodal[b_index] >= 2)
                b_outcome = 1;

            if (gshare[g_index] >= 2)
                g_outcome = 1;

            if (hybrid[h_index]<2)
            {
                if(actual_outcome == 1 && bimodal[b_index] < 3)
                    bimodal[b_index]++;
                else if (actual_outcome == 0 && bimodal[b_index] > 0)
                    bimodal[b_index]--;

                if (b_outcome != actual_outcome)
                    mispredictions++;
            }
            else
            {
                if (actual_outcome == 1 && gshare[g_index] < 3)
                    gshare[g_index]++;
                else if (actual_outcome == 0 && gshare[g_index] > 0)
                    gshare[g_index]--;

                if (g_outcome != actual_outcome)
                    mispredictions++;
            }

            gbh = gbh>>1;
            int temp = actual_outcome<<(N-1);
            gbh = gbh | temp;

            if (g_outcome == actual_outcome)
            {
                if (b_outcome != actual_outcome)
                {
                    if (hybrid[h_index] < 3)
                        hybrid[h_index]++;
                }
            }
            else
            {
                if (b_outcome == actual_outcome)
                {
                    if (hybrid[h_index] > 0)
                        hybrid[h_index]--;
                }
            }
        }
    }
    //Calculate the index value of the given address
    void calculate_index(unsigned long int addr)
    {
        int gbh_mask, index_mask, hybrid_mask, temp_mask;
        int difference = M1 - N;
        temp_mask = (pow(2,difference) - 1);

        if(mode == 0)     
        {   //Bimodal
            index_mask = pow(2,M2) - 1;
            b_index = addr & index_mask;
        }
        else if(mode == 1)     
        {   //Gshare
            index_mask = pow(2,M1) - 1;
            gbh_mask = pow(2,N) - 1;
            g_index = addr & index_mask;
            
            int temp = g_index & temp_mask;

            g_index = (g_index>>difference)^(gbh & gbh_mask);
            g_index = (g_index<<difference)|temp;
        }
        else if(mode == 2)     
        {   //Hybrid
            index_mask = pow(2,M2) - 1;
            b_index = addr & index_mask;

            index_mask = pow(2,M1) - 1;
            gbh_mask = pow(2,N) - 1;       
            g_index = addr & index_mask;

            int temp = g_index & temp_mask;

            g_index = (g_index>>difference)^(gbh & gbh_mask);
            g_index = (g_index<<difference)|temp;

            hybrid_mask = pow(2,K) - 1;
            h_index = addr & hybrid_mask;
        }
    }
    //Print prediction table contents
    void print_prediction_table()
    {
        if(mode == 0)
        {
            cout<<"\nFINAL BIMODAL CONTENTS ";
            int x = pow(2,M2);
            for (int i = 0; i < x; i++)
                cout<<"\n "<<i<<" "<<bimodal[i];
        }
        else if (mode == 1)
        {
            cout<<"\nFINAL GSHARE CONTENTS ";
            int x = pow(2,M1);
            for (int i = 0; i < x; i++)
                cout<<"\n "<<i<<" "<<gshare[i];
        }
        else if (mode == 2)
        {
            cout<<"\nFINAL CHOOSER CONTENTS ";
            int x = pow(2,K);
            for (int i = 0; i < x; i++)
                cout<<"\n "<<i<<" "<<hybrid[i];

            cout<<"\nFINAL GSHARE CONTENTS ";
            int y = pow(2,M1);
            for (int i = 0; i < y; i++)
                cout<<"\n "<<i<<" "<<gshare[i];

            cout<<"\nFINAL BIMODAL CONTENTS ";
            int z = pow(2,M2);
            for (int i = 0; i < z; i++)
                cout<<"\n "<<i<<" "<<bimodal[i];
        }
    }
};
#endif
