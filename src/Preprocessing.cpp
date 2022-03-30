
#include "stdc++.h"
#include <vector>
#include <boost/dynamic_bitset.hpp> 
#include <bitset>
#include <map>
#include <fstream>
using namespace std; 
using namespace boost;

char complement(char c)
{
    if (c == 'A') return 'T';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    if (c == 'T') return 'A';
    return c;
}
int validNucleotide( vector<char> seq){
    std::string str (seq.begin(), seq.end());
    std::size_t found = str.find_first_not_of("ACTG");
    if (found!=std::string::npos)
        return 1;
    return 0;
}
static void appendLineToFile(string filepath, string line)
{
    std::ofstream file;
    file.open(filepath, std::ios::out | std::ios::app);
    if (file.fail())
        throw std::ios_base::failure(std::strerror(errno));

    //make sure write fails with exception if something is wrong
    file.exceptions(file.exceptions() | std::ios::failbit | std::ifstream::badbit);

    file << line << std::endl;
}
vector<char> reversecomplement(vector<char> seq){
    vector<char> temp;
    temp.clear();
    for (auto const& ntide : seq){
        temp.push_back(complement(ntide));
    }
    std::reverse(temp.begin(),temp.end());
    return temp;
}

template <size_t N1, size_t N2 >
bitset <N1 + N2> concat( const bitset <N1> & b1, const bitset <N2> & b2 ) {
    string s1 = b1.to_string();
    string s2 = b2.to_string();
    return bitset <N1 + N2>( s1 + s2 );
}
boost::dynamic_bitset<uint8_t> twoBitsEncoding( vector <char> seq ){
    //std::vector<int> ohe;
    //ohe.clear();
    std::reverse(seq.begin(), seq.end());
    boost::dynamic_bitset<uint8_t> ohe;
    //string ohe_string = "";
    
    for (auto const& ntide : seq){
        
        
        if (ntide == 'A'){
                //ohe_string += "11";
                ohe.push_back(1);
                ohe.push_back(1);
                ohe.push_back(1);
                ohe.push_back(1);

            }
        else if (ntide == 'C')
                {
                //ohe_string += "01";
                ohe.push_back(1);
                ohe.push_back(1) ;
                ohe.push_back(0);
                ohe.push_back(1) ; 
                }
        else if (ntide == 'T')
                {
                //ohe_string += "00";
                ohe.push_back(0);
                ohe.push_back(0);
                ohe.push_back(0);
                ohe.push_back(1) ; 
                }
         else if (ntide == 'G')
                {
                //ohe_string += "10";
                ohe.push_back(1);
                ohe.push_back(0) ; 
                ohe.push_back(1);
                ohe.push_back(0) ; 
                }
            //else
              //  {
             //   ohe_string += "000";
                
              //  }
    }
        
    //dynamic_bitset<uint8_t> ohe(ohe_string);
    //std::string seq_string (seq.begin(), seq.end());
    //cout << ohe << "   " << seq_string << endl;
    return ohe;
}
int SizeTToInt(size_t data)
{
    if (data > std::numeric_limits<int>::max())
        throw std::runtime_error("Invalid cast.");
    return static_cast<int>(data);
}
std::vector<string> readToBuffer(char* buffer, size_t size) {
    //int count_of_reads = 0;
    int count = 0;
    std::vector<size_t> posa;
    std::vector<string> seqs;
    

    //get the number of lines and end of line position
    const char * p = &buffer[0];

    for(size_t is=0; is<size; is++) {
        if(p[is] == '\n' ) { //should add EOF as well to be safe || buffer[i]!= EOF
            count++;
            posa.push_back(is);
        }
    }

    //#pragma omp parallel for
    //count = 500;
    for(int i=1; i<count ;i++) { 
        const int len = SizeTToInt(posa[i] - posa[i-1]);
        char* buff = &buffer[posa[i-1]];
        std::string line(buff, len);
        std::size_t found = line.find('>');

        if (found==std::string::npos){
            
            seqs.push_back(line);
        }
        //count_of_reads++;
        //printf("\n Total No. of reads: %d , %d \n",count_of_reads,  count);


           
    }

    printf("\n Finished reading the file! \n");
    posa.clear();
    return seqs;
}

boost::dynamic_bitset<uint8_t> patternEncoding11( vector <char> seq ){

    std::vector<char> A_string ;
    std::vector<char> C_string ;
    std::vector<char> T_string ;
    std::vector<char> G_string ;
    
    for (auto const& ntide : seq){
        
        
        if (ntide == 'A'){
                A_string.push_back('1');
                C_string.push_back('0');
                T_string.push_back('0');
                G_string.push_back('0');

            }
        else if (ntide == 'C')
                {
                A_string.push_back('0');
                C_string.push_back('1');
                T_string.push_back('0');
                G_string.push_back('0');
                
                }
        else if (ntide == 'T')
                {
                A_string.push_back('0');
                C_string.push_back('0');
                T_string.push_back('1');
                G_string.push_back('0');
                
                }
         else if (ntide == 'G')
                {
                A_string.push_back('0');
                C_string.push_back('0');
                T_string.push_back('0');
                G_string.push_back('1');
                 
                }
            else
                {
                A_string.push_back('0');
                C_string.push_back('0');
                T_string.push_back('0');
                G_string.push_back('0');
                
                }
    }
    std::vector<char> phash (A_string.size()* 4);
    
    std::move(A_string.begin(), A_string.end(), phash.begin());
    std::move(C_string.begin(), C_string.end(), phash.begin()+ A_string.size());
    std::move(T_string.begin(), T_string.end(), phash.begin()+ A_string.size()*2);
    std::move(G_string.begin(), G_string.end(), phash.begin()+ A_string.size()*3);
    boost::dynamic_bitset<uint8_t> phash_bitset;
    for (int i=phash.size()-1; i>=0; i--){
        phash_bitset.push_back(phash[i]-'0');
    }
    //cout<<phash_bitset.size()<<" "<< sizeof(phash_bitset) <<endl;
    return phash_bitset;
}
boost::dynamic_bitset<uint8_t> OHEncoding( vector <char> seq ){

    std::vector<char> bit_string ;

    
    for (auto const& ntide : seq){
        
        
        if (ntide == 'A'){
                bit_string.push_back('1');
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');

            }
        else if (ntide == 'C')
                {
                bit_string.push_back('0');
                bit_string.push_back('1');
                bit_string.push_back('0');
                bit_string.push_back('0');
                
                }
        else if (ntide == 'T')
                {
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('1');
                bit_string.push_back('0');
                
                }
         else if (ntide == 'G')
                {
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('1');
                 
                }
            else
                {
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');
                
                }
    }

    boost::dynamic_bitset<uint8_t> vecbit;
    for (int i=bit_string.size()-1; i>=0; i--){
        vecbit.push_back(bit_string[i]-'0');
    }
    //cout<<vecbit.size()<<" "<< sizeof(vecbit) <<endl;
    return vecbit;
}
boost::dynamic_bitset<uint8_t> patternEncoding( vector <char> seq ){

    std::vector<char> bit_string ;

    
    for (auto const& ntide : seq){
        
        
        if (ntide == 'A'){
                bit_string.push_back('1');
                bit_string.push_back('1');
                bit_string.push_back('1');
                bit_string.push_back('1');

            }
        else if (ntide == 'C')
                {
                bit_string.push_back('1');
                bit_string.push_back('1');
                bit_string.push_back('1');
                bit_string.push_back('0');
                
                }
        else if (ntide == 'T')
                {
                bit_string.push_back('1');
                bit_string.push_back('1');
                bit_string.push_back('0');
                bit_string.push_back('0');
                
                }
         else if (ntide == 'G')
                {
                bit_string.push_back('1');
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');
                 
                }
            else
                {
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');
                bit_string.push_back('0');
                
                }
    }

    boost::dynamic_bitset<uint8_t> vecbit;
    for (int i=bit_string.size()-1; i>=0; i--){
        vecbit.push_back(bit_string[i]-'0');
    }
    //cout<<vecbit.size()<<" "<< sizeof(vecbit) <<endl;
    return vecbit;
}
void vecVec2Arr(std::vector<std::vector<uint8_t>> vectors , uint8_t *arrcode ){
    int total_size = vectors.size()* vectors[0].size();

    // 2. Create a vector to hold the data.
    std::vector<uint8_t> flattened;
    flattened.reserve(total_size);

    // 3. Fill it
    for (auto& vec : vectors)
        for (std::uint8_t elem : vec){
            flattened.push_back(elem);
            //cout<< "vector filled with" << int(elem)<< endl;
        }
    cout << total_size<< endl;
    cout << flattened.size()<< endl;
    for (int i=0; i<20; i++){
        cout << (int)flattened[i] << std::endl;
    } 
    // 4. Obtain the array
    //std::copy(flattened.begin(), flattened.end(), arr);
    arrcode = flattened.data();
    cout << sizeof(arrcode)/sizeof(*arrcode)<< endl;
    for (int i=0; i<20; i++){
        cout << (int)arrcode[i] << std::endl;
    } 
    
}
std::vector<uint8_t> splituint8t( const dynamic_bitset<uint8_t> bits )
{
    
    std::vector< uint8_t > result ;
    result.reserve(bits.num_blocks());
    boost::to_block_range( bits, std::back_inserter(result) ) ;
    return result;
}
std::map <int, std::vector<std::vector<uint8_t>> > createRefs(int n , string filepath ){
    std::map <int, std::vector<std::vector<uint8_t>> > refMap;
    std::vector<std::vector<uint8_t>> parsrefvec;
    string homedir = getenv("HOME");
    string revfile = homedir + "/Data/batcovidKReads/reads_(" + filepath + ")_1.fq";
    string frwfile = homedir + "/Data/batcovidKReads/reads_(" + filepath + ")_2.fq";
    ifstream input1(revfile.c_str());
    ifstream input2(frwfile.c_str());
    if(!input1.good() || !input2.good()){
        cerr << "Error opening " << revfile <<". Bailing out." << endl;
        //return {-1};
    }
    string line1, name1, content1,line2, name2, content2;
    int linnum = 1;
    int cnt = 0;
    while( getline( input1, line1 ).good() && getline( input2, line2 ).good()){
        if (linnum % 4 == 2){
            cnt++;
            string frag = line1+line2;
            vector <char> Genchars(frag.begin(), frag.end());
            
            for (int i = 1; i <= (Genchars.size() - n) ; ++i) {
                
                vector <char> subList1 (Genchars.begin() + i, Genchars.begin() + i + n );
                if (validNucleotide(subList1)==0){
                    vector <char> revcom1 = reversecomplement(subList1);

                    const dynamic_bitset<uint8_t> ohe_sublist1  = OHEncoding(subList1);

                    std::vector<uint8_t> sublist_bytes1;
                    sublist_bytes1 = splituint8t(ohe_sublist1);
                    parsrefvec.push_back(sublist_bytes1);

                    const dynamic_bitset<uint8_t> ohe_revcom1  = OHEncoding(revcom1);
                    
                    std::vector<uint8_t> revcom_bytes1;
                    revcom_bytes1 = splituint8t(ohe_revcom1);
                    parsrefvec.push_back(revcom_bytes1);
                    subList1.clear();
                    revcom1.clear();
                }
            }

            refMap[cnt] = parsrefvec;
            

            linnum ++;
            frag.clear();
            parsrefvec.clear();
        }
        else linnum ++;
    }
    input1.close();
    input2.close();
    cout<<refMap.size() << "| " <<refMap[1].size()<<" "<< sizeof(refMap[1]) <<endl;

    return refMap;
}
std::vector<std::vector<uint8_t>> createRefs_spike(int n , string filepath ){
    std::vector<std::vector<uint8_t>> parsrefvec;
    string seqs="";
    
    string homedir = getenv("HOME");
    string revfile = homedir + "/Data/" + filepath + ".fa";
    ifstream input1(revfile.c_str());
    if(!input1.good()){
        cerr << "Error opening " << revfile <<". Bailing out." << endl;
        //return {-1};
    }
    string line1, name1, content1;
    while( ! input1.eof()){
        getline( input1, line1 );
        std::size_t found = line1.find('>');

        if (found==std::string::npos){
            seqs+=line1;
        }
        
    }

    vector <char> Genchars(seqs.begin(), seqs.end());

    for (int i = 1; i <= (Genchars.size() - n) ; ++i) {

        vector <char> subList1 (Genchars.begin() + i, Genchars.begin() + i + n );
        if (validNucleotide(subList1)==0){
            vector <char> revcom1 = reversecomplement(subList1);

            const dynamic_bitset<uint8_t> ohe_sublist1  = OHEncoding(subList1);

            std::vector<uint8_t> sublist_bytes1;
            sublist_bytes1 = splituint8t(ohe_sublist1);
            parsrefvec.push_back(sublist_bytes1);

            const dynamic_bitset<uint8_t> ohe_revcom1  = OHEncoding(revcom1);

            std::vector<uint8_t> revcom_bytes1;
            revcom_bytes1 = splituint8t(ohe_revcom1);
            parsrefvec.push_back(revcom_bytes1);
            subList1.clear();
            revcom1.clear();
        }

    }
    input1.close();
    return parsrefvec;
}
std::vector<std::vector<uint8_t>> createdb(int n , std::vector<string> batpoops ){
    std::set<std::vector<uint8_t>> parsset;
    for (auto file: batpoops)
    {
        string f = "/home/pjsheini/Data/LargeDB-SRAs/"+ file +".fa";
        
        const char* filepath = f.c_str();
        printf("\n %s \n",filepath);
        cout<< "strat creating db" << endl;
        // read db fasta file
        FILE * pFile;
        long long lSize;
        char * buffer;
        size_t fileSize;
        

        pFile = fopen ( filepath, "rb" );
        if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

        // obtain file size:
        fseek (pFile , 0 , SEEK_END);
        lSize = ftell (pFile);
        rewind (pFile);
        
        // allocate memory to contain the whole file:
        buffer = (char*) malloc ((sizeof(char)*lSize));
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

        // copy the file into the buffer:
        fileSize = fread (buffer,1,lSize,pFile);
        if (fileSize != lSize) {fputs ("Reading error",stderr); exit (3);}

        /* the whole file is now loaded in the memory buffer. */
        std::vector<string> seqreads = readToBuffer(buffer, fileSize);
        // terminate
        fclose (pFile);
        free (buffer);
        //string p2f = "/home/pjsheini/C-Jupyter/hwt2_OptaneV01/kmers_db.txt";
        //#pragma omp parallel for //shared (parsset) 
        for( int iseq=0; iseq<seqreads.size(); iseq++)
        {
            vector <char> Genchars(seqreads[iseq].begin(), seqreads[iseq].end());
            #pragma omp parallel for shared (parsset)
            for (int i = 1; i <= (Genchars.size() - n) ; ++i) {
                vector<char> subList (Genchars.begin() + i, Genchars.begin() + i + n );
                //std::string kmstring (Genchars.begin() + i, Genchars.begin() + i + n );
                //appendLineToFile(p2f,kmstring);
                if (validNucleotide(subList)==0){
                    vector <char> revcom = reversecomplement(subList);

                    const dynamic_bitset<uint8_t> ohe_sublist  = OHEncoding(subList);
                    std::vector<uint8_t> sublist_bytes;
                    sublist_bytes = splituint8t(ohe_sublist);
                    #pragma omp critical
                    parsset.insert(sublist_bytes);

                    const dynamic_bitset<uint8_t> ohe_revcom  = OHEncoding(revcom);
                    std::vector<uint8_t> revcom_bytes;
                    revcom_bytes = splituint8t(ohe_revcom);
                    #pragma omp critical
                    parsset.insert(revcom_bytes);

                    sublist_bytes.clear();
                    revcom_bytes.clear();
                    subList.clear();
                    revcom.clear();
                }

            }

        }
    }
    std::vector<std::vector<uint8_t>> parsvec(parsset.begin(), parsset.end());

    cout << parsvec.size()<< "vectorSize"<<endl;
    cout<<parsvec[0].size()<<" "<< typeid(parsvec[0]).name() <<endl;
    return parsvec;
}
