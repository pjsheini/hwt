#include "root_node.h"
#include "preprocessing.h"
#include "io.h"
#include "nodeops.h"
#include "bitops.h"
#include "SearchNode.h"
#include "memusage.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "node.h"
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
int Node::B_over_8;
uint8_t* Node::dbcode;
struct timespec begin1, end1, begin0, end0, begin11, end11, begin00, end00;
timespec diff(timespec start, timespec end )
{
    timespec temp;

    if ((end.tv_nsec-start.tv_nsec)<0)
    {
            temp.tv_sec = (end.tv_sec-start.tv_sec-1);
            temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    }
    else 
    {
            temp.tv_sec = (end.tv_sec-start.tv_sec);
            temp.tv_nsec = (end.tv_nsec-start.tv_nsec);
    }
    return temp;
}
//std::vector<char> Node::dbcode;
int main(){
    std::vector<string> batFeces = {"BatCov_RpYN06.unmapped", "BatCov_RmYN05.unmapped","BatCov_fb_HlYN18_unmapped"};
    //std::vector<string> batFeces = {"COVID_19"};
    string dbpooppath = "/home/pjsheini/Data/bat.unmapped.fa";
    string dbpath  = "/home/pjsheini/Data/COVID_19.fa";
    int n=32;
    //int R = 3;
    int B = 128;
    UINT32 K =10;  //the number of nearest neighbors to be retrieved
    //UINT32 Q =0 ; //the number of query points to use from <infile>, default all
                //N = Set the number of binary codes from the beginning of the dataset file to be used
                //-T <number>          Set the maximum size of each node in the tree
    string Query_Genome = "RaTG13" ; //#Ebola, "RaTG13","RmYN01"
    int B_over_8 = B/8;
    int capacity;
    uint8_t *codes_db;
    int dim1codes;
    uint8_t *codes_query;
    int dim1queries;
    printf("Loading codes... ");
    fflush(stdout);
    //std::array<uint8_t,(size_t)N * (B/8)> codes_db;
    std::vector<std::vector<uint8_t>> parsvec = createdb( n, batFeces);
    UINT32 N = parsvec.size();
    codes_db = (uint8_t*)malloc((size_t)N * (B/8) * sizeof(uint8_t));
    //cout<< (size_t)N <<endl;
    ////////////////////
    //vecVec2Arr(parsvec, codes);
    int total_size = parsvec.size()* parsvec[0].size();

    // 2. Create a vector to hold the data.
    std::vector<uint8_t> flattened;
    flattened.reserve(total_size);

    // 3. Fill it
    for (auto& vec : parsvec)
        for (std::uint8_t elem : vec){
            flattened.push_back(elem);
        }

    codes_db = flattened.data();
    flattened.clear(); 
    ////////////////////////
    printf("Loaded\n");
    fflush(stdout);



    capacity = 5;
    Node* curr_node;
    uint8_t* ccode;
    int tree_lvl = 0;
    int max_tree_lvl = 0;

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin11);
    clock_gettime(CLOCK_REALTIME, &begin00);

    Node* root_node = new Node(-1,NULL);
    root_node->isleaf = false;
    Node::B_over_8 = B_over_8;
    Node::dbcode = codes_db;
    int max_lvl = 10;
    int expansions=0;
    printf("Started creating tree\n");
    fflush(stdout);
    for(uint64_t i=0;i<N;i++)
    {

        tree_lvl = 0;


        ccode = codes_db + i*B_over_8;
        //cout<< (int)ccode[1] << (int)ccode[2]<<endl;
        curr_node = root_node;
        while(true){
            if(curr_node->isleaf) {
                    
                curr_node->push_to_node(ccode,i);
                //cout<< curr_node << " | ";
                if(curr_node->size == capacity && tree_lvl<max_lvl){
                    curr_node->expand();
                    expansions++;
                }
                break;
            }
            else{
                uint8_t subnorms[pow2(tree_lvl)];
                norm_chunks(subnorms,tree_lvl,ccode,B_over_8);
                curr_node = curr_node->find_the_child(ccode,subnorms);
            }

            tree_lvl++;
        }
        if(tree_lvl>max_tree_lvl)
            max_tree_lvl = tree_lvl;
        if(i> N-1)
        {
            printf("%ld items are loaded\n",i);
            fflush(stdout);
        }


    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end11);
    clock_gettime(CLOCK_REALTIME, &end00);
    auto wt = diff(begin00, end00 );
    auto cput = diff( begin11,end11);
    printf("cput = %ld seconds and %ld nanoseconds, wall = %ld seconds and %ld nanoseconds\n",cput.tv_sec , cput.tv_nsec, wt.tv_sec, wt.tv_nsec);
    printf("Done Loading\n");

    double vm = -1;
    double rss = -1;
    process_mem_usage(&vm, &rss);
        vm  /= double(1024);
        rss /= double(1024);
        printf("VM %.1fmb | RSS %.1fmb\n",vm,rss);
    printf("max tree level = %d\n",max_tree_lvl);
    printf("expansions = %d\n",expansions);
    printf("Done\n");
    std::map<int, std::vector<std::vector<uint8_t>> > refvec = createRefs( n, Query_Genome);
    for (int R=1; R<11;R++){
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin1);
        clock_gettime(CLOCK_REALTIME, &begin0);
        UINT32 NQ = 10;//refvec[1].size();
        codes_query = (uint8_t*)malloc((size_t)NQ * (B/8) * sizeof(uint8_t));
        uint8_t *codesq, *ccodeq;
        printf("Loading %d queries\n",NQ);

        //////////////////////////
        //vecVec2Arr(parsrefvec, codesq);
        //int totalref_size = parsrefvec.size()* parsrefvec[0].size();
        int totalref_size = refvec.size()* refvec[1][2].size();
        int fragfound = 0 ;
        int sumkmfound = 0;
        int allvisitednodes = 0;
        uint64_t total_num_comp=0;
        for (auto const& [fragnum, parsrefvec]: refvec){
            // 2. Create a vector to hold the data.
            std::vector<uint8_t> flattenedrefs;
            flattenedrefs.reserve(totalref_size);
            if (fragnum>1)
                break;
            // 3. Fill it
            for (auto& vec : parsrefvec)
                for (std::uint8_t elem : vec){
                    flattenedrefs.push_back(elem);
                }

            codes_query = flattenedrefs.data();
            ccodeq = codesq = codes_query;
            sumkmfound = 0;
            allvisitednodes = 0;
            total_num_comp=0;
            flattenedrefs.clear();
            SearchNode* sn = new SearchNode(codes_db, B_over_8,max_tree_lvl,R,K, root_node,NQ);


            for (UINT32 i=0; i < NQ; i++) {
                sn->setQuery(ccodeq);

                ccodeq += B_over_8;

                auto [found , nodevisit] = sn->HNN_search();
                if (found>0)
                    sumkmfound += 1;

                allvisitednodes+=nodevisit;
                total_num_comp += sn->num_compare;


            }
            delete sn;
        }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end1);
        clock_gettime(CLOCK_REALTIME, &end0);
        auto wttt = diff(begin0, end0 );

        auto cputt = diff( begin1,end1);
        

        printf("number of comparisons = %ld\n",total_num_comp);
        printf("Search Finished");
        cout<< Query_Genome;
        printf(" %d  %d outOf ( %d )   %ld.%ld \n",R-1, sumkmfound,  NQ,cputt.tv_sec , cputt.tv_nsec);
    }
    delete root_node;
}
