//
//  main.cpp
//  memoryEfficientReferenceGenomeFMIndex
//
//  Created by muhammed oguzhan kulekci on 8.10.2019.
//  Copyright © 2019 muhammed oguzhan kulekci. All rights reserved.
//

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/lcp.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include <fstream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <iomanip>
#include <divsufsort.h>
#include <chrono>
#include <random>

#include <pthread.h>

//default sdsl values are SA_SAMPLE_FREQ=32 ISA_SAMPLE_FREQ=64
#define SA_SAMPLE_FREQ  4
#define ISA_SAMPLE_FREQ 64
#define MYDEBUGx
#define BLCDx

using namespace sdsl;
using namespace std;

//#define SA_SAMPLING_STRAT sa_order_sa_sampling<>
#define SA_SAMPLING_STRAT text_order_sa_sampling<>

#ifdef BLCD
typedef csa_wt<wt_blcd<rrr_vector<127>>,SA_SAMPLE_FREQ,ISA_SAMPLE_FREQ,SA_SAMPLING_STRAT, isa_sampling<>> fmindexType;
#else
//typedef csa_wt<wt_hutu<bit_vector_il<>,rank_support_il<>>,SA_SAMPLE_FREQ,ISA_SAMPLE_FREQ,SA_SAMPLING_STRAT, isa_sampling<>> fmindexType;
typedef csa_wt<wt_hutu<bit_vector,rank_support_v5<>>,SA_SAMPLE_FREQ,ISA_SAMPLE_FREQ,SA_SAMPLING_STRAT, isa_sampling<>> fmindexType;
#endif

typedef struct cacheItem{
    unsigned long long left;
    unsigned long long right;
    unsigned long long count;
} cacheItemType;

typedef struct threadParameter{
    unsigned long long begin;
    unsigned long long end;
    unsigned long long visitedSymbolCount;
    unsigned long long scannedSymbolCount;
} threadParameterType;


fmindexType fm_index;
lcp_dac<> lcp;

unsigned char* finalBWT;
char *buffer;
unsigned int  MAXCACHEITEMS;
pthread_t* threads;
unsigned char complement[256];

unsigned long long fileSize(char* fname){
    std::ifstream in(fname, std::ifstream::ate | std::ifstream::binary);
    unsigned long int filesize = in.tellg();
    in.close();
    return filesize;
}

unsigned long int prepareSequence(const char* fileName,unsigned long long filesize){
    std::ifstream in(fileName, std::ifstream::in | std::ifstream::binary);
    in.read(&buffer[0],filesize);
    in.close();
    cout << "\nFile Name and File Size in bytes\t:" << fileName << " " << filesize << endl;
    unsigned char* readPtr  = (unsigned char*)&buffer[0]; //readPtr is guaranteed to be larger than or equal to the writePtr
    unsigned char* writePtr = (unsigned char*)&buffer[0];
    while (readPtr<(unsigned char*)&buffer[filesize]){
        switch (*readPtr){
            case 'A':
            case 'a': *writePtr='A'; writePtr++; break;
            case 'C':
            case 'c': *writePtr='C'; writePtr++; break;
            case 'G':
            case 'g': *writePtr='G'; writePtr++; break;
            case 'T':
            case 't': *writePtr='T'; writePtr++; break;
            case 'N':
            case 'n': *writePtr='E'; writePtr++; break; //replace Ns with Es to keep up the line in the complementing as A<C<E<G<T and T>G>E>C>A
        }
        readPtr++;
    }
    *writePtr = '\0';  writePtr++;
    return writePtr - (unsigned char*)&buffer[0]; //length of the DNA sequence including the NULL at the end thus GATTACA\0 is a 8 symbols sequence
}

long int  visit(unsigned long long *left, unsigned long long *right, unsigned char symbol, fmindexType* fmindex, unsigned long long endOfSeqPosition){ //HEAVY DUTY
    long int count=0;
    unsigned char complementedSymbol = complement[symbol];
    auto res = fmindex->wavelet_tree.lex_count(*left, (*right)+1, complementedSymbol);
    count = get<2>(res);// number of symbols that are larger than the symbol
    if ((*left <= endOfSeqPosition) && (*right >= endOfSeqPosition)) count++; // if end of sequence marker sign is in between the boundaries increase the lexicographically smaller suffix count.
    backward_search(*fmindex, *left, *right, complementedSymbol, *left, *right); // update the interval
    return count;
}

void  insertSymbol(unsigned char symbol, unsigned long long position, unsigned char* T){
    T[position] = symbol;
    if ('\0'==symbol) T[position] = '$';
}

void* doWork(void* tparam){
    
    unsigned long long startPosition = ((threadParameterType*)tparam)->begin;
    unsigned long long endPosition = ((threadParameterType*)tparam)->end;
    
    cacheItemType* cache = new cacheItemType[MAXCACHEITEMS+1];
    //fm_index = cst->csa;
    
    unsigned char symbol;
    unsigned long long p=0,C=0,rankReverse=0,L=0,R = fm_index.size()-1,dollarSignPos = fm_index.isa[0],saValue;
    
    unsigned long long visitedSymbols = 0;
    
    for(unsigned long long s= startPosition;s>endPosition;s--){
        
        //printProgress((double)(startPosition-s) / (double)(startPosition-endPosition) );
        
        visitedSymbols = 0;
        saValue = fm_index[s];
        
        while((long)L<=(long)R){
            
            if (p<=MAXCACHEITEMS){ //save current range and rank
                cache[p].left=L;
                cache[p].right=R;
                cache[p].count=rankReverse;
            }
            
            symbol = buffer[saValue+p];
            if (symbol=='\0'){
                if ((L<=dollarSignPos) && (R>=dollarSignPos)){
                    rankReverse++; // this suffix on forward strand is also a suffix on the reverse complement strand. As end of sequence is smaller than the forward-reverse splitter, add 1 to the rank
                }
                break;
            }else {
                C = visit(&L,&R,symbol,&fm_index,dollarSignPos);
                rankReverse += C;
                visitedSymbols++;
            }
            p++;
        }// while
#ifdef MYDEBUG
        cout << s << ' ' << rankReverse << ' ' << rankReverse+s << ' ' << fm_index.bwt[s] << ' ' << visitedSymbols << ' ' << lcp[s] << endl;
#endif
        insertSymbol(fm_index.bwt[s],rankReverse + s,finalBWT);
        //if (visitedSymbols>10000)
        //cout << s << ' ' << fm_index[s] << ' ' << visitedSymbols << ' ' << lcp[s] << endl;
        
        if ( p <= lcp[s]){
            
            do{ //no need to look for the previous suffixes while their lcp values are greater then p since their initial strings of length p+1 are always unique
                s--;
                saValue = fm_index[s];
                insertSymbol(buffer[saValue-1],rankReverse + s,finalBWT); //fm_index.bwt[s]
                //visitedSymbols = 0;
                //cout << s << ' ' << visitedSymbols << ' ' << cst.lcp[s] << endl;
            }while ((p<=lcp[s]) && (s>endPosition));
        }
        
        p = lcp[s]; // guaranteed that p>cst.lcp[s]
        if (p>MAXCACHEITEMS) p = MAXCACHEITEMS;
        
        L = cache[p].left;
        R = cache[p].right;
        rankReverse = cache[p].count;
    }
    cout << visitedSymbols << endl;
    ((threadParameterType*)tparam)->visitedSymbolCount = visitedSymbols;
    return tparam;
}


void* doWorkEfficient(void* tparam){
    
    unsigned long long startPosition = ((threadParameterType*)tparam)->begin;
    unsigned long long endPosition   = ((threadParameterType*)tparam)->end;
    
    cacheItemType* cache = new cacheItemType[MAXCACHEITEMS+1];
    
    unsigned char symbol;
    unsigned long long offset=0, C=0, rankReverse=0, L=0, R = fm_index.size()-1, dollarSignPos = fm_index.isa[0],saValue;
    bool stop;
    
    unsigned long long visitedSymbols = 0;
    unsigned long long scannedSymbols = 0;
    
    cache[0].left  = L;
    cache[0].right = R;
    cache[0].count = C;
    
    offset = 0;
    
    for(unsigned long long s= startPosition;s>endPosition;s--){ //for each suffix
        
        saValue = fm_index[s]; //starting position of this suffix on the text
        
        stop = false;
        
        L = cache[offset].left;
        R = cache[offset].right;
        rankReverse = cache[offset].count;
        
        while(!stop){
            if (L>R){
                stop=true;
            }else if (L<R){ // perform HEAVY DUTY on the position
                
                visitedSymbols++;
                
                symbol = buffer[saValue+offset];
                
                if (symbol=='\0'){
                    
                    if ((L<=dollarSignPos) && (R>=dollarSignPos)){
                        rankReverse++; // if this happens, it means this suffix on forward strand is also a suffix on the reverse complement strand. As end of sequence is smaller than the forward-reverse splitter, add 1 to the rank
                    }
                    stop = true; // no deed to go further after observing the smallest symbol, notice no need to update the cache since we didnot perform and lexCount or backwards search
                    
                }else {
                    C = visit(&L,&R,symbol,&fm_index,dollarSignPos); // count smaller suffixes on reverse complement
                    rankReverse += C; // accumulate this to the rankReverse
                    offset++; //increase the offset
                    //if (offset<MAXCACHEITEMS){ // if we have enough cache size reserved
                        cache[offset].left  = L; //save the values in the cache for future use
                        cache[offset].right = R;
                        cache[offset].count = rankReverse;
                    //}
                }
                
            }else{// means L==R, then perform light scanning
                unsigned char * LeftToRightPtr = (unsigned char*) &buffer[saValue+offset];
                unsigned char * RightToLeftPtr = (unsigned char*) &buffer[fm_index[L]-1];
                unsigned char complementedSymbol;
                while(1){
                    scannedSymbols++;
                    if (RightToLeftPtr < (unsigned char*)&buffer[0]){
                        rankReverse++;
                        stop=true;
                        break;
                    }else if (LeftToRightPtr > (unsigned char*)&buffer[fm_index.size()-2]){
                        stop=true;
                        break;
                    }else{
                        complementedSymbol = complement[*RightToLeftPtr];
                        if (*LeftToRightPtr<complementedSymbol){
                            stop = true;
                            break;
                        } else if (*LeftToRightPtr > complementedSymbol){
                            stop = true;
                            rankReverse++;
                            break;
                        }
                    }
                    LeftToRightPtr++;
                    RightToLeftPtr--;
                }
            }//else L=R
        }// while
        
        insertSymbol(fm_index.bwt[s],rankReverse + s,finalBWT);
        if (offset>lcp[s]){
            offset = lcp[s];
        }
        if (offset>MAXCACHEITEMS) offset = MAXCACHEITEMS;
#ifdef MYDEBUG
        cout << s << " final Rank: " << rankReverse + s << endl;
#endif

    }
    
    ((threadParameterType*)tparam)->visitedSymbolCount = visitedSymbols;
    ((threadParameterType*)tparam)->scannedSymbolCount = scannedSymbols;

    return tparam;
}


int main(int argc, const char * argv[]) {
    //Usage: mergi input_file_name number_of_threads
    
    complement['A'] = 'T'; complement['C'] = 'G'; complement['E'] = 'E'; complement['G'] = 'C'; complement['T'] = 'A';
    unsigned int numberofthreads = atoi(argv[2]);
    threads = new pthread_t[numberofthreads];
    std::chrono::high_resolution_clock::time_point t1,t2;
    
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // read file and load into memory...
    t1 = std::chrono::high_resolution_clock::now();
    unsigned long long fSize =  fileSize((char*)argv[1]);
    buffer = new char[fSize+1];
    memset(buffer, 0, fSize+1);
    unsigned long int seqSize = prepareSequence(argv[1],fSize); //seqSize is the number of bases in the forward reference genome including the NULL symbol at the end
    t2 = std::chrono::high_resolution_clock::now();
    cout << "File preperation in seconds:\t"<<   std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // allocate space for the finalBWT result
    finalBWT = new unsigned char[2*seqSize+1];
    memset(finalBWT,0,2*seqSize+1); //reset to NULL
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // construct the fm_index on forward text
    t1 = std::chrono::high_resolution_clock::now();
    cache_config cc(false);
    construct_im(fm_index, (const char*)buffer,1);
    t2 = std::chrono::high_resolution_clock::now();
    cout << "T_1_^1: Forward strand FM-index construction time in seconds:\t"<<   std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // compute the LCP array for efficient calculation
    t1 = std::chrono::high_resolution_clock::now();
    cc.delete_files=true;
    construct_im(lcp,(const char*)buffer, 1);
    t2 = std::chrono::high_resolution_clock::now();
    cout << "T_1^2: Forward strand LCP array computation time in seconds:\t"<<  std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    
    cout << "S_2: FM-index + LCP size in MB: " << size_in_mega_bytes(fm_index) + size_in_mega_bytes(lcp) << " = " << size_in_mega_bytes(fm_index) << " + " << size_in_mega_bytes(lcp) << endl;

    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // find the maximum length LCP
    t1 = std::chrono::high_resolution_clock::now();
    unsigned long long maxLCP=0;
    for(unsigned long long h=0;h<lcp.size();h++) if (lcp[h]>maxLCP) maxLCP = lcp[h];
    MAXCACHEITEMS = maxLCP;
    t2 = std::chrono::high_resolution_clock::now();
    cout << "Finding maximum LCP value (= " << maxLCP << ") time in seconds:\t"<<  std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    cout << "Cache memory in MB " << numberofthreads * (maxLCP * 24)/(1024*1024) << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
#ifdef MYDEBUG
    cout << " i SA ISA  BWT  T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3B %:3T", fm_index);
    cout << lcp << endl;
#endif

    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // prepare the threads to run doWork
    t1 = std::chrono::high_resolution_clock::now();
    threadParameterType parameter[numberofthreads];
    threadParameterType *out[numberofthreads];
    unsigned long long blockLength = (seqSize-1) / numberofthreads;
    parameter[0].begin = seqSize-1;
    parameter[0].end   = seqSize-blockLength;
    for(unsigned int t=1; t<numberofthreads;t++){
        parameter[t].begin  = parameter[t-1].end;
        parameter[t].end    = parameter[t-1].end - blockLength;
    }
    parameter[numberofthreads-1].end = 0;
    for(unsigned int t=0; t<numberofthreads;t++){
        pthread_create(&threads[t], NULL, doWorkEfficient, (void*)&parameter[t]);
        //cout << endl << parameter[t].begin << '\t' << parameter[t].end << " STARTED."  << endl;
    }
    for(unsigned int t=0; t<numberofthreads;t++){//wait threads to finish their duties
        pthread_join(threads[t], (void**)(&out[t]));
        cout << "RETURNED: " << out[t]->begin << ' ' << out[t]->end << ' ' << out[t]->visitedSymbolCount << ' ' << out[t]->scannedSymbolCount << endl;
    }
    finalBWT[1] = buffer[seqSize-2]; // the symbol jsut preceding the inter sequence marker #, i.e. the last symbol of F in F#R$, notice that we assume  $<#<A<C<E(N)<G<T
    t2 = std::chrono::high_resolution_clock::now();
    cout << "T_2: Rank calculations of forward suffixes time in seconds:\t"<<   std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //reverse strand = reverse complement of the input sequence
    t1 = std::chrono::high_resolution_clock::now();
    unsigned long long middle = (seqSize%2==0)?seqSize/2:(seqSize-1)/2;
    for(unsigned long long i=0; i< middle ; i++){
        unsigned char temp;
        temp = complement[(unsigned char)buffer[i]];
        buffer[i] = complement[buffer[seqSize-2-i]];
        buffer[seqSize-2-i] = temp;
    }
    finalBWT[0] = buffer[seqSize-2];
    t2 = std::chrono::high_resolution_clock::now();
    cout << "Reverse complementing the forward strand in seconds:\t"<<   std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Reverse strand FM-Index computation
    t1 = std::chrono::high_resolution_clock::now();
    construct_im(fm_index, (const char*)buffer,1);
    t2 = std::chrono::high_resolution_clock::now();
    cout << "T_3: Reverse strand FM-index construction time in seconds:\t"<<  std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl ;
    cout << "S_3: Reverse strand FM-index size:" << size_in_mega_bytes(fm_index) << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
#ifdef MYDEBUG
    cout << " i SA ISA   BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S  %2p  %3B  %:3T", fm_index);
#endif

    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //Locating the reverse strand suffixes into the final BWT
    t1 = std::chrono::high_resolution_clock::now();
    unsigned long long posfinalBWT=0;
    while (finalBWT[posfinalBWT]!='\0') posfinalBWT++;
    for(unsigned long int s=1;s<fm_index.size();s++){
        finalBWT[posfinalBWT] = fm_index.bwt[s];
        if ('\0'==finalBWT[posfinalBWT]) finalBWT[posfinalBWT] = '#';
        posfinalBWT++;
        while (finalBWT[posfinalBWT]!='\0') posfinalBWT++;
    }
    t2 = std::chrono::high_resolution_clock::now();
    cout << "Reverse strand suffix placement into final BWT string in seconds:\t"<<  std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl ;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    

    /* VERIFICATION PROCESS: check whether the computed finalBWT is equal to the BWT constructed over the whole reference genome */
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // compute the whole forwardstrand#reversestrand FM-Index
    delete[] buffer;
    t1 = std::chrono::high_resolution_clock::now();
    fSize =  fileSize((char*)argv[1]);
    buffer = new char[2*seqSize+1];
    memset(buffer, 0, fSize+1);
    prepareSequence(argv[1],fSize); //read file into buffer
    buffer[seqSize-1]='#';
    char* tmpptr=&buffer[seqSize];
    for(long long i=seqSize-2; i>=0 ; i--){
        *tmpptr = complement[buffer[i]];
        tmpptr++;
    }
    *tmpptr='\0';
    t2 = std::chrono::high_resolution_clock::now();
    cout << "VERIFICATION:\nRead file and append the reverse complementing time in seconds:\t"<<   std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl;
    
    t1 = std::chrono::high_resolution_clock::now();
    construct_im(fm_index, (const char*)buffer,1);
    t2 = std::chrono::high_resolution_clock::now();
    cout << "T_BASE: Forwar+Reverse strand FM-index construction time in milliseconds:\t"<<  std::chrono::duration_cast<std::chrono::seconds>( t2 - t1 ).count() << endl ;
    
#ifdef MYDEBUG
    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", fm_index);
#endif
    
    cout << "Verifying the correctness of the generated BWT ...." << endl;
    for(unsigned long int s=2;s<fm_index.size();s++)
        if ((fm_index.bwt[s]!=finalBWT[s]) && (fm_index.bwt[s]!='\0')) cout << "ERROR at position " << s << ' ' << fm_index.bwt[s] << ' ' << finalBWT[s] << endl;
    cout << "DONE" << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
#ifdef MYDEBUG
    cout << finalBWT << endl;
#endif
    return 0;
}
