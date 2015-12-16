/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Jason Warta
 * 2015-12-15
 * CS 301: Assembly Language Programming
 * 
 * Huffman Compression Algorithm 
 * 
 * Due to ineffeiciencies in the buffer system I used, this program has known
 * issues compressing and decompressing files over 100 KB or files with abnormally
 * high numbers of some characters
 * 
 * to adjust the size of the header depending on the number of characters
 * in the file, change the following typedef:
 * freq_type : type for counting character frequency, currently a uint32_t
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
#include <fstream>
using std::fstream;
using std::ios;
using std::ios_base;
#include <sstream>
using std::stringstream;
#include <string>
using std::string;
#include <algorithm>
using std::sort_heap;
using std::sort;
using std::swap;
#include <vector>
using std::vector;
#include <deque>
using std::deque;
#include <cstdint>
using std::uint64_t;
using std::uint32_t;
using std::uint16_t;
using std::uint8_t;
using std::size_t;
#include <map>
using std::map;
#include <cmath>
using std::max;
using std::min;
#include <utility>
using std::pair;
#include <thread>
using std::thread;
#include <queue>
using std::queue;
#include <functional>
using std::ref;

//using bitset rather than defining my own class for bit handling
// because manually manipulating bits proved to be somewhat error-prone
#include <bitset> 
using std::bitset;

//adjust the type of freq_type and totBits_type to increase or decrease filesize
typedef uint32_t freq_type; // for counting character frequency
typedef uint8_t tail_type; // for tracking number of extra bits at end of file
typedef uint64_t totBits_type; // for tracking total bits in a file
typedef uint8_t magic_type; //type of magic number

//max chars, used in various loops
const int CHARS = 256;

//fingerprint to determine if a file was encoded by this program
const magic_type MAGIC = 'c' + 's';

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * H_Node struct
 * used by Huffman compression algorithm when mapping a file
 * keeps track of two H_Node pointers, the character, and associated frequency
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct H_Node{
	H_Node(char data, freq_type freq=0):data_(data),freq_(freq){}
	H_Node(H_Node * left, H_Node * right, freq_type freq=0):left_(left),right_(right),freq_(freq){}

	H_Node * left_ = nullptr;
	H_Node * right_ = nullptr;
	unsigned char data_;
	freq_type freq_ = 0;
};//end struct H_Node

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * PointerCompare()
 * allows comparison of two H_Node pointers
 * pass to std::sort as the comparsion function
 * 	when sorting a container of H_Node pointers
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct PointerCompare{
	bool operator() (const H_Node * lhs, const H_Node * rhs){
		return (*lhs).freq_ < (*rhs).freq_;
	}
};//end PointerCompare struct

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * operator<
 * overloaded operator for H_Node objects
 * evaluates equality of H_Node nodes based on their frequency value
 * if two H_Node nodes have equal frequency, the left node is marked greater
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool operator<( H_Node & lhs, H_Node & rhs ) { 
	if(lhs.freq_ == rhs.freq_){
		if(lhs.left_ == nullptr && lhs.right_ == nullptr){
			return true;
		}

		if(rhs.left_ == nullptr && rhs.right_ == nullptr){
			return false;
		}
	}
	return lhs.freq_ < rhs.freq_; 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * RANGE struct, instantiated as 'range'
 * track min and max code lengths
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
struct RANGE{
  uint8_t min=8;
  uint8_t max=1;
} range;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * BytesToFreq
 * union for converting a block of bytes to character frequency
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
typedef union{
	char bytes[sizeof(freq_type)];
	freq_type num;
} BytesToFreq;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * invertNum(BytesToFreq &)
 * swaps the bytes of a block of chars to convert to a big-Endian number
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void invertNum(BytesToFreq & btf){
	swap(btf.bytes[0],btf.bytes[3]);
	swap(btf.bytes[1],btf.bytes[2]);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * printEncoding(ofstream &, map<unsigned char,string>)
 * writes the map portion of the header to the ofstream
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void printEncoding(fstream & ofs, map<unsigned char,string> & codeMap){
	char map_size[sizeof(uint8_t)];
	*map_size = codeMap.size();
	ofs.write(map_size,sizeof(uint8_t));

	char character[1] = {};
	char length[1] = {};
	char code[sizeof(freq_type)] = {};

	for(auto i : codeMap){
		*character = i.first;
		ofs.write(character,1);

		*length = i.second.length();
		ofs.write(length,1);

		bitset<64> bits(i.second);
		*code = (freq_type)(bits.to_ullong());
		ofs.write(code,sizeof(freq_type));
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * readEncoding(ifstream &, map<unsigned char, string)
 * reads the map portion of the header from the ifstream
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void readEncoding(fstream & ifs, map<unsigned char, string> & codeMap){
	//get the number of nodes in the map
	char size[sizeof(uint8_t)];
	ifs.read(size,sizeof(uint8_t));
	uint8_t map_size = *size;

	//set up buffers
	char character[1] = {};
	char length[1] = {};
	char code[sizeof(freq_type)] = {};

	unsigned char d_char = 0;
	unsigned long tempLong = 0;
	uint8_t len = 0;
	string temp_str;

	for(size_t i = 0; i < map_size; i++){
		//get the character
		ifs.read(character,1);
		d_char = character[0];

		//get the length of the code
		ifs.read(length,1);
		len = *length;

		//update the stored value for min and max code lengths
		range.min = min( len, range.min );
		range.max = max( len, range.max );

		//get the code
		ifs.read(code,sizeof(freq_type));
		tempLong = (*code);
		bitset<64> code = tempLong;
		temp_str = code.to_string();
		//trim code, amd map the char and accompanying code
		codeMap[d_char] = temp_str.substr(temp_str.size()-len,temp_str.size());
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * countChars(freq_type *, ifstream &)
 * Counts the frequency of characters in a file
 * iterates through the ifstream and tallies the characters
 * the value in the array at index "char" is incremented for each character found
 * 
 * used by compression algorithm
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void countChars(freq_type * charCount, fstream & fs){
	unsigned char * buff = new unsigned char[1];

	fs.seekg(ios_base::beg);
	while(fs.good()){
		fs.read( (char*)buff,1); //read one byte at a time
		charCount[ *buff ]++;
	}
	buff = NULL;
	delete buff;
	
	fs.seekg(ios_base::beg);
	fs.clear();
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * printCharCount(freq_type *, ofstream &)
 * iterates through an array of character frequencies and prints them to the file
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void printCharCount(freq_type * charCount, fstream & fs){
	BytesToFreq buff;
	
	for(int i = 0; i < CHARS; i++){
		buff.num = charCount[i];
		invertNum(buff);
		fs.write( buff.bytes, sizeof(freq_type) );
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * printMap(map<,> &)
 * templated function, prints out an std::map
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
template<typename A, typename B>
void printMap(map<A,B> & input){
	for(auto i : input){
		cout << i.first << " " << i.second << endl;
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * readMagic(ifstream &)
 * reads the magic number from the file
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool readMagic(fstream & fs){
	uint64_t temp = 0;
	char buff[8];
	fs.read(buff,sizeof(uint64_t));
	for(int i = sizeof(uint64_t); i >= 0; i--){
		temp = (temp << 8) + buff[i];
	}
	return temp == MAGIC;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * readCharCount(freq_type *, ifstream &)
 * reads the character frequency from a compressed file
 * 
 * used by decompression
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void readCharCount(freq_type * charCount, fstream & fs){
	BytesToFreq buff;

	for(int i = 0; i < CHARS; i++){
		fs.read( buff.bytes, sizeof(freq_type) );
		invertNum(buff);
		charCount[i] = buff.num;
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * buildTree(H_node * &, freq_type *)
 * takes int ptr to array with tallies of char counts
 * builds a binary tree based on the frequency of the characters in the file
 * 
 * used by both comression and decompression
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void buildTree(H_Node * & tree, freq_type * charCount){
	deque<H_Node*> freqList;

	//file deque with nodes for each character
	for(int i = 0; i < CHARS ; i++){
		if( charCount[i] ){
			// if the frequency of a character is greater than 0, add it to the vector
			freqList.push_back(new H_Node( (unsigned char)i, charCount[i] ) );
		}
	}

	//sort the deque
	sort(freqList.begin(), freqList.end(), PointerCompare() );

	while(freqList.size() > 1){
		
		//create a new node point that is the root of the first two items in the deque, 
		//and has the combined frequency of those two items
		freqList.push_back(new H_Node( freqList[0], freqList[1], freqList[0]->freq_ + freqList[1]->freq_ ) );
		
		//remove the two elements just used
		freqList.pop_front();
		freqList.pop_front();

		//resort the deque
		sort(freqList.begin(), freqList.end(), PointerCompare() );
	}

	tree = freqList[0];
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * getBit(unsigned char, int)
 * takes a byte and a position
 * gets the bit at the given position in the char
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
size_t getBit(unsigned char byte, int position){
	//get bit at given position and return
	return (byte >> position) & 0x1;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * setBit(char &, int, size_t)
 * takes a b
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void setBit(char & byte, int position, size_t val){
	//code copied from stackoverflow
	if(getBit(byte,position) == 0){
		byte |= 1 << position;
	}
	else if(getBit(byte,position) == 1){
		byte &= ~(1 << position);
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * mapCodes(H_Node *, map<unsigned char, string> &, string)
 * takes a H_Node pointer to the root of a tree
 * takes a reference to a map of unsigned chars and strings
 * takes a string. this should be emtpy. defaults to blank
 * maps the tree by recursive transversal
 * each time a leaf is found, the char and its path (code) are added to the map
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void mapCodes(H_Node * tree, map<unsigned char,string> & codeMap, string path=""){
	if(tree->left_ == nullptr && tree->right_ == nullptr){
		codeMap[tree->data_] = path;
	} else {
		mapCodes(tree->left_,codeMap,path+"0");
		mapCodes(tree->right_,codeMap,path+"1");
	}
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * printEncodedByte(ofstream &, deque<size_T> &)
 * prints a block of 8 'bits' from the passed deque
 * items are removed from the queue after being printed
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void printByte(fstream & fs, deque<size_t> & bitStream){
	unsigned char byte;

	//put the first bit in the byte
	byte + bitStream[0];
	//delete first bit from queue
	bitStream.pop_front();

	for(int i = 0; i < 7 && i <= bitStream.size(); i++){
		//delete bit from queue
		bitStream.pop_front();
	}
	fs << byte;
}

/*//currently non-functional, originally written for multithreaded version of compression
void readPlainText(fstream & ifs, map<unsigned char, string> & codeMap, queue<size_t> & bitStream, bool & complete){
	if(ifs.is_open()){

	}
	// ifs.open("file.txt",ios::in | ios::binary);
	// cout << "arrived in input function" << endl;
	unsigned char byte = 0;
	// cout << "checking on file" << endl;
	ifs.seekg(ios_base::beg);
	if(ifs.is_open()){
		// cout << "file is open" << endl;
		ifs.clear();
		while(!ifs.eof()){
			ifs.read( (char*)( (long)byte), 1 );
			// cout << byte << endl;

			string code = codeMap[byte];
			// cout << "doing stuff 1" << endl;
			for(size_t i = 0; i < code.size(); i++){
				bitStream.push( code[i]-48 );
				// cout << "doing stuff 2" << endl;
			}
		}
	} else {
		cout << "file i/o error when reading plain text" << endl;
	}
	
	complete = true;
}
//currently no-functional, originally written for mulithreaded version of decompression
void writeEncodedBytes(fstream & ofs, queue<size_t> & bitStream, bool & complete){
	// cout << "arrived in output function" << endl;
	unsigned char byte = 0;
	// byte[0] = 0;
	size_t bit = 8;

	if(ofs.is_open()){
		cout << "file is open" << endl;
		cout << complete << endl;
		while(!complete || (complete && bitStream.size() > 0)){
			cout << "processing bits. stream size: " << bitStream.size() << endl;
			if(bitStream.size() > 0){
				cout << "setting bit " << bit << endl;
				// setBit(byte,bit-1,bitStream.front());
				bitStream.pop();
				bit--;
				cout << "bit " << bit << " has been written" << endl;
			}

			if(bit == 0){
				ofs.write( (char*)( (long)byte), 1 );
				bit = 8;
			}
		}

		if(bit != 8){
			cout << "appending trailing bits" << endl;
			tail_type extraBits = bits;
			// extraBits = bit;

			while(bit != 0){
				// setBit(byte,bit-1,0);
				bit--;
			}
			ofs << byte;
			// ofs.write( (char*)byte,1 );

			ofs.seekp(ios_base::beg);
			// ofs << extraBits;
			ofs.write( (char*)( (long)extraBits),1 );
			cout << "done with trailing bits" << endl;
		}
	} else {
		cout << "file i/o error when writing encoded bytes" << endl;
	}
}*/

/* //currently non-functioning, written originally for multithreaded version of decompression
void readEncodedBytes(fstream & ifs, map<string,unsigned char> & charMap, queue<unsigned char> & byteStream,bool & complete){
	unsigned char byte[1];
	stringstream code;
	if(ifs.is_open()){
		while(ifs.good()){
			ifs.read( (char*)byte, 1 );

			for(int i = 7; i >= 0; i++){
				code << getBit(*byte,i);
			}

			for(int i = 0; i < code.str().size(); i++){
				if( charMap.find( code.str().substr(0,i) ) != charMap.end() ){
					byteStream.push(charMap[code.str().substr(0,i)]);
					code.str().erase(0,i);
					i = 0;
				}
			}

		}		
	}
	complete = true;
}

//currently non-functional, written originally for multithreaded version of decompression
void writePlainText(fstream & ofs, queue<unsigned char> & byteStream,bool & complete){
	if(ofs.is_open()){
		while(!complete || (complete && byteStream.size() > 0)){
			if(byteStream.size() > 0){
				ofs << byteStream.front();
				// ofs.write( (char*)byteStream.front(), 1 );
				byteStream.pop();
			}
		}
	}
}*/

/*
 * encodeChars
 * reads file a char at a time, adds char code to queue by calling getCode
 * prints char by calling printEncodedByte
 * 
 * used by compression algorithm
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * encodeChars(ifstream &, ofstream &, map<unsigned char,string> &)
 * takes input and output filestreams
 * reads input filestream and prints encoded bits to the output filestream
 * returns true if succesfully completed
 * returns false if there was an error with either filestream
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool encodeFile(fstream & ifs, fstream & ofs, map<unsigned char, string> & codeMap){
	if(ifs.is_open() && ofs.is_open()){

		//write fingerprint to start of file
		char magic[sizeof(magic_type)];
		bitset<64> bits(MAGIC);
		*magic = (magic_type)(bits.to_ullong());
		ofs.write(magic,sizeof(magic_type));

		//set placeholder for extraBits value
		char placeholder[sizeof(totBits_type)];
		*placeholder = 0;
		ofs.write(placeholder,1);

		printEncoding(ofs,codeMap);

		//set file position and clear any bad flags
		ifs.seekg(ios_base::beg);
		ifs.clear();

		//set variables for later use
		char iBuff[1];
		string codeBuff = "";
		char extraBits[1];
		size_t bytes = 0;

		//read contents of file to string
		while(ifs.good()){
			ifs.read(iBuff,1);
			codeBuff += codeMap[*iBuff];
			// extraBits += codeMap[*iBuff].length();
		}

		//add additional bits to code buffer to bring the last byte up to 8 bits
		if(codeBuff.length()%8 != 0){
			extraBits[0] = 8 - (codeBuff.length() % 8);
			for(uint8_t i = 0; i < extraBits[0]; i++){
				codeBuff += "0";
			}
		}

		bytes = codeBuff.length()/8;
		char oBuff[1];

		//work through the code buffer, converting string of 8 'bits' to bytes and writing them to the file
		for(int i = 0; i < codeBuff.length()-1; i+=8){
			bitset<8> bits(codeBuff.substr(i,i+8));
			oBuff[0] = (uint8_t)(bits.to_ulong());
			ofs.write(oBuff,1);
		}

		// *placeholder = totalBits
		//jump back to beginning of file and write the number of extra bits
		ofs.seekp(ios_base::beg + sizeof(MAGIC));
		ofs.write( extraBits, 1 );
		// ofs.write( placeholder, sizeof(totBits_type) );
		return true;
	} else {
		cout << "File IO error in encodeFile function." << endl;
		return false;
	}
	return true;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * decodeFile(ifstream &, ofstream &, map<string, unsigned char> &, tail_type)
 * provide input and output filestreams
 * provide a map with strings containing codes as the key, unsigned char as the mapped value
 * provide a number type with the number of extra bits at the end of the file
 * reads encoded file from input filestream, writes decoded bytes to output filestream
 * returns true if completed successfully
 * returns false if there were errors with the filestream
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool decodeFile(fstream & ifs, fstream & ofs, map<string, unsigned char> & charMap, tail_type trailingBits){
	if(ifs.is_open() && ofs.is_open()){
    char iBuff[1];
    char oBuff[1];
    string sBuff;
    string val;

    while(ifs.good()){
    	ifs.read(iBuff,1);
      bitset<8> bits(iBuff[0]);
      sBuff += bits.to_string();

      if(sBuff.length() > range.min && ifs.peek() != EOF)

      if(ifs.peek() == EOF) break;
    }

    sBuff.erase( sBuff.length() - trailingBits, sBuff.length() );

    for(uint64_t i = range.min; i < sBuff.length(); i++){
    	val = sBuff.substr(0,i);
    	if(charMap.find(val) != charMap.end()){
    		oBuff[0] = charMap[val];
    		ofs.write(oBuff,1);
    		sBuff.erase(0,val.length());
    		i=0;
    	}
    }
  }else {
  	cout << "File IO error in decodeFile function" << endl;
  	return false;
  }
  return true;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * compression(string,string)
 * provide input and output filenames
 * compresses the input file and stores it in the output file
 * returns true if successfully completed
 * returns false if there was an error
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool compression(string inFile, string outFile){
	fstream ifs(inFile, ios::in | ios::binary | ios::ate);
	fstream ofs(outFile, ios::out | ios::binary | ios::ate);

	if(ifs.is_open() && ofs.is_open()){
		ifs.seekg(ios_base::beg);
		ofs.seekp(ios_base::beg);

		//count char frequency in file
		freq_type charCount[CHARS] = {};
		countChars(charCount,ifs);

		//build tree from char freqency list
		H_Node * tree;
		buildTree(tree,charCount);

		map<unsigned char, string> codeMap;
		mapCodes(tree,codeMap,"");

		bool complete = false;
		queue<size_t> bitStream;

		bool success = encodeFile(ifs,ofs,codeMap);

		// thread output(writeEncodedBytes,ref(ofs),ref(bitStream),ref(complete));
		// readPlainText(ifs,codeMap,bitStream,complete);
		// output.join();

		ifs.close();
		ofs.close();
		return success;
	} else {
		return false;
	}
	return true;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * decompression(string,string)
 * provide input and output filenames
 * decompresses the input file and stores it in the output file
 * returns true if successfully completed
 * returns false if there was an error
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool decompression(string inFile, string outFile){
	fstream ifs(inFile, ios::in | ios::binary | ios::ate);
	fstream ofs(outFile, ios::out | ios::binary | ios::ate);

	if(ifs.is_open() && ofs.is_open()){
		ifs.seekg(ios_base::beg);
		ofs.seekp(ios_base::beg);

		//get magic number
		magic_type magicNum;
		char magic[sizeof(magic_type)];
		ifs.read(magic,sizeof(magic_type));
		magicNum = (magic_type)(*magic);
		//check magic number, terminate function if it doesn't match
		if(magicNum != MAGIC){
		 	cout << "That file was not encded by this program." << endl;
		 	return false;
		}

		//get number of extra bits
		char buff[1];
		tail_type extraBits;
		ifs.read( buff, sizeof(tail_type) );
		extraBits = (tail_type)(*buff);

		//get map from comressed file
		map<unsigned char, string> codeMap;
		readEncoding(ifs,codeMap);

		//reverse the map's indexing
		map<string,unsigned char> charMap;
		for(auto iter = codeMap.begin(); iter != codeMap.end(); iter++){
			charMap[iter->second] = iter->first;
		}

		bool complete = false;
		queue<unsigned char> byteStream;

		bool success = decodeFile(ifs,ofs,charMap,extraBits);

		// thread output(writePlainText,ref(ofs),ref(byteStream),ref(complete));
		// readEncodedBytes(ifs,charMap,byteStream,complete);
		// output.join();

		ifs.close();
		ofs.close();
		return success;
	} else {
		return false;
	}
	return false;
}
	
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * mapFile(string)
 * given a file name, creates and maps the huffman tree of a file
 * prints out the characters with their accompaning code
 * 
 * this is for demonstration purposes only
 * it does not perform any compression or decompressions, and is
 * not a required step to compress or decompress a file
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
bool mapFile(string inFile){
	fstream ifs(inFile, ios::in | ios::binary | ios::ate);

	if(ifs.is_open()){
		ifs.seekg(ios_base::beg);

		//get char frequency from compressed file
		freq_type charCount[CHARS] = {};
		countChars(charCount,ifs);

		//build tree from char frequency list
		H_Node * tree;
		buildTree(tree,charCount);

		//map the tree
		map<unsigned char, string> codeMap;
		mapCodes(tree,codeMap,"");

		printMap(codeMap);
		return true;
	}
	else {
		return false;
	}
	return false;
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * inputError()
 * prints instructions for program usage to std::cout
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void inputError(){
	cout << "\n";
	cout << "Incorrect usage.\n";
	cout << "Run the program using the following options.\n";
	cout << "  Compression:\n";
	cout << "   -c SOURCE DESTINATION\n";
	cout << "  Decompression:\n";
	cout << "   -d SOURCE DESTINATION\n";
	cout << "  Map file: only for deomnstration purposes\n";
	cout << "   -m SOURCE\n";
	cout << endl;
}

int main(int argc, char** argv){
	//arg 0 is the command used to launch the app, ie "./a.out" or "valgrind ./a.out"
	//subsequent args are the options

	if(argc == 4){
		//get args
		string option(argv[1]);
		string ifname(argv[2]);
		string ofname(argv[3]);

		if(option == "-c"){
			if(compression(ifname,ofname)) cout << "Finished compressing file." << endl;
			else cout << "Encountered error when compressing." << endl;
		} 
		else if(option == "-d"){
			if(decompression(ifname,ofname)) cout << "Finished decompressing file." << endl;
			else cout << "Encountered error when decomressing." << endl;
		} 
		else {
			inputError();
		}
	} else if(argc == 3){
		//get args
		string option(argv[1]);
		string ifname(argv[2]);

		if(option == "-m"){
			if(mapFile(ifname)) cout << "Finished mapping file." << endl;
			else cout << "Encountered error when mapping file." << endl;
		} 
		else {
			inputError();
		}
	} 
	else {
		inputError();
	}

	return 0;
}