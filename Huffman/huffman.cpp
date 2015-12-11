
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
using std::uint32_t;//for tracking character count
using std::uint64_t;//magic number
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
#include <bitset>
using std::bitset;
// #include "boost/lockfree/queue.hpp"

typedef uint64_t freq_type;
typedef uint8_t tail_type;

typedef union{
	char bytes[sizeof(freq_type)];
	freq_type num;
} BytesToFreq;

//struct for tracking range of code lengths
struct RANGE{
  uint8_t min=8;
  uint8_t max=1;
} range;

const int CHARS = 256;
const uint64_t MAGIC = 861314319301;

unsigned char bits;
size_t bitCount = 0;

struct H_tree{
	H_tree(char data, freq_type freq=0):data_(data),freq_(freq){}
	H_tree(H_tree * left, H_tree * right, freq_type freq=0):left_(left),right_(right),freq_(freq){}
	H_tree(H_tree * left, H_tree * right, unsigned char data):left_(left),right_(right),data_(data){}

	H_tree * left_ = nullptr;
	H_tree * right_ = nullptr;
	unsigned char data_;
	freq_type freq_ = 0;
};//end struct tree

/*
 * PointerCompare function
 * compares two pointers to H_tree for sorting
 */
struct PointerCompare{
	bool operator() (const H_tree * lhs, const H_tree * rhs){
		return (*lhs).freq_ < (*rhs).freq_;
	}
};

/*
 * operator<
 * overloaded operator for H_tree objects
 * evaluates equlity of H_tree nodes based on their frequency value
 * if two H_tree nodes have equal frequency, 
 */
bool operator<( H_tree & lhs, H_tree & rhs ) { 
	// cout << "arrived in operator<" << endl;
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

void invertNum(BytesToFreq & btf){
	swap(btf.bytes[1],btf.bytes[4]);
	swap(btf.bytes[2],btf.bytes[3]);
}

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

void readEncoding(fstream & ifs, map<unsigned char, string> & codeMap){
	char size[sizeof(uint8_t)];
	ifs.read(size,sizeof(uint8_t));
	uint8_t map_size = *size;

	char character[1] = {};
	char length[1] = {};
	char code[sizeof(freq_type)] = {};

	unsigned char d_char = 0;
	unsigned long tempLong = 0;
	uint8_t len = 0;
	string temp_str;

	for(size_t i = 0; i < map_size; i++){
		ifs.read(character,1);
		d_char = character[0];

		ifs.read(length,1);
		len = *length;

		range.min = min( len, range.min );
		range.max = max( len, range.max );

		ifs.read(code,sizeof(freq_type));
		tempLong = (*code);
		bitset<64> code = tempLong;
		temp_str = code.to_string();
		codeMap[d_char] = temp_str.substr(temp_str.size()-len,temp_str.size());
	}
}

/* * * * * * * * *
 *  Compression  *
 * * * * * * * * */

/*
 * countChars
 * takes int pointer initialized to const CHARS items
 * takes fstream by reference
 * tallys the number of characters by incrementing that index of the int array
 */
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

/*
 * printMagic
 * takes writable fstream by reference
 * prints the magic number to the fstream
 */
void printMagic(fstream & fs){
	fs << MAGIC;
}

/*
 * printCharCount
 * takes array of freq_type containing char frequencies
 * takes writable fstream by reference
 * prints char frequencies to the fstream
 */
void printCharCount(freq_type * charCount, fstream & fs){
	BytesToFreq buff;
	
	for(int i = 0; i < CHARS; i++){
		buff.num = charCount[i];
		invertNum(buff);
		fs.write( buff.bytes, sizeof(freq_type) );
	}
	// char buff[ sizeof(freq_type) ];
	// for(int i = 0; i < CHARS; i++){
	// 	*buff = charCount[i];
	// 	fs.write( buff, sizeof(freq_type) );
	// }
}

template<typename A, typename B>
void printMap(map<A,B> input){
	for(auto i : input){
		cout << i.first << " " << i.second << endl;
	}
}

/*
 * printTrailingBits
 * takes freq_type containing total number of bits needed to compress the file
 * takes writeable fstream by reference
 * prints out the number of leftover bits to the file
 */
// void printTrailingBits(freq_type freq, fstream & fs){
// 	fs << (tail_type)(8 - (freq % 8) );
// }

/* * * * * * * * * *
 *  Decompression  *
 * * * * * * * * * */

/*
 * readMagic
 * read number from file, return true if it is the magic number
 */
bool readMagic(fstream & fs){
	uint64_t temp = 0;
	char buff[8];
	fs.read(buff,sizeof(uint64_t));
	for(int i = sizeof(uint64_t); i >= 0; i--){
		temp = (temp << 8) + buff[i];
	}
	return temp == MAGIC;
}

/*
 * readCharCount
 * reads the character frequency array from an encoded file
 */
void readCharCount(freq_type * charCount, fstream & fs){
	BytesToFreq buff;

	for(int i = 0; i < CHARS; i++){
		fs.read( buff.bytes, sizeof(freq_type) );
		invertNum(buff);
		charCount[i] = buff.num;
	}
}

/*
 * readTrailingBits
 * reads the number indicating the number of trash bits at the end of the file
 */
tail_type readTrailingBits(fstream & fs){
	tail_type tail;
	fs.read( (char*)((long)tail),sizeof(tail_type) );
	return tail;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * buildTree
 * takes int ptr to array with tallies of char counts
 * builds a binary tree based on the frequency of the characters in the file
 * 
 * used by both comression and decompression
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void buildTree(H_tree * & tree, freq_type * charCount){
	deque<H_tree*> freqList;

	//file deque with nodes for each character
	for(int i = 0; i < CHARS ; i++){
		if( charCount[i] ){
			// if the frequency of a character is greater than 0, add it to the vector
			freqList.push_back(new H_tree( (unsigned char)i, charCount[i] ) );
		}
	}

	sort(freqList.begin(), freqList.end(), PointerCompare() );

	while(freqList.size() > 1){
		
		//create a new node point that is the root of the first two items in the deque, 
		//and has the combined frequency of those two items
		freqList.push_back(new H_tree( freqList[0], freqList[1], freqList[0]->freq_ + freqList[1]->freq_ ) );
		//remove the two elements just used
		freqList.pop_front();
		freqList.pop_front();

		sort(freqList.begin(), freqList.end(), PointerCompare() );
	}

	tree = freqList[0];
}

size_t getBit(unsigned char byte, int position){
	//get bit at given position and return
	return (byte >> position) & 0x1;
}

void setBit(char & byte, int position, size_t val){
	//code copied from stackoverflow
	if(getBit(byte,position) == 0){
		byte |= 1 << position;
	}
	else if(getBit(byte,position) == 1){
		byte &= ~(1 << position);
	}
}

/*
 * mapCodes
 * takes a H_tree pointer
 * takes a map<unsigned char, string> by reference
 * takes a string
 * maps the tree by recursive transversal
 * each time a leaf is found, the char and path to it are added to the map
 */
void mapCodes(H_tree * tree, map<unsigned char,string> & codeMap, string path){
	if(tree->left_ == nullptr && tree->right_ == nullptr){
		codeMap[tree->data_] = path;
	} else {
		mapCodes(tree->left_,codeMap,path+"0");
		mapCodes(tree->right_,codeMap,path+"1");
	}
}

/*
 * printEncodedByte
 * build a char using the first 8 bits in the queue
 * print the created char
 * 
 * used by compression algorithm
 */
void printByte(fstream & fs, deque<size_t> & bitStream){
	unsigned char byte;

	//put the first bit in the byte
	byte + bitStream[0];
	//delete first bit from queue
	bitStream.pop_front();

	for(int i = 0; i < 7 && i <= bitStream.size(); i++){

		// setBit(byte,7-i,bitStream[0]);
		//left shift the byte by 1
		// byte = byte << 1;
		//add next bit to the byte
		// byte + bitStream[0];
		//delete bit from queue
		bitStream.pop_front();
	}
	fs << byte;
}

/*
 * getChars
 * reads from ifs 
 */
// unsigned char getChars(fstream & ifs, fstream & ofs, H_tree * tree){
// 	unsigned char byte;

// 	deque<size_t> bitStream;

// 	if(ifs.good()){

// 		while(ifs.good() && !(tree->left_ == nullptr && tree->right_ == nullptr) ){
// 			ifs.read( (char)byte, 1 );

// 			for(int i = 7; i >= 0; i--){
// 				bitStream.push_back( getBit(byte,i) );
// 			}

// 			if(bitStream.size() >= 8){
// 				printEncodedByte(ofs,bitStream);
// 			}
// 		}
// 	}
// }
void readPlainText(fstream & ifs, map<unsigned char, string> & codeMap, queue<size_t> & bitStream, bool & complete){
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
}

/*
 * encodeChars
 * reads file a char at a time, adds char code to queue by calling getCode
 * prints char by calling printEncodedByte
 * 
 * used by compression algorithm
 */
bool encodeFile(fstream & ifs, fstream & ofs, map<unsigned char, string> & codeMap){
	if(ifs.is_open() && ofs.is_open()){

		//set placeholder for 
		char placeholder[1];
		placeholder[0] = 0;
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
			// if(ifs.peek() == EOF) break;
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

		//jump back to beginning of file and write the number of extra bits
		ofs.seekp(ios_base::beg);
		ofs.write( extraBits, 1 );

	} else {
		cout << "File IO error in encodeFile function." << endl;
		return false;
	}
	return true;
}

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
}

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

      if(ifs.peek() == EOF) break;
    }

    sBuff.erase( sBuff.length() - trailingBits, sBuff.length() );

    for(int i = 0; i < sBuff.length(); i++){
    	val = sBuff.substr(0,i);
    	if(charMap.find(val) != charMap.end()){
    		oBuff[0] = charMap[val];
    		ofs.write(oBuff,1);
    		sBuff.erase(0,val.length());
    		i=0;
    	}
    }
  }
/*

    while(ifs.good()){
      ifs.read(iBuff,1);
      encodedText.push_back(iBuff);
      
      bitset<8> byte(*iBuff);
      for(int i = 0; i < 8; i++){
        sBuff += byte.to_string();
      }
      
      eof = (ifs.peek() == EOF);
      
      //if we are at EOF
      if(eof){
        sBuff.erase(sBuff.length()-trailingBits,sBuff.length());
      }
      
      if(sBuff.length() >= range.min){
        for(int i = range.min; i <= range.max && i < sBuff.length(); i++){
          if( charMap.find( sBuff.substr(0,i) ) != charMap.end() ){
            obuff[0] = charMap[sBuff.substr(0,i)];
            sBuff.erase(0,i);
            ofs.write(obuff,1);
          }
        }
      }
    }
  } else {
  	return false;
  }	*/
  return true;
}

void printStuff(string stuff){
	cout << stuff << endl;
}

/*
 * compression
 */
bool compression(string inFile, string outFile){
	fstream ifs(inFile, ios::in | ios::binary | ios::ate);
	fstream ofs(outFile, ios::out | ios::binary | ios::ate);

	if(ifs.is_open() && ofs.is_open()){
		ifs.seekg(ios_base::beg);
		ofs.seekp(ios_base::beg);

		// char placeholder[1];
		// placeholder[0] = 0;
		// ofs.write( placeholder, sizeof(tail_type) );

		//count char frequency in file
		freq_type charCount[CHARS] = {};
		countChars(charCount,ifs);

		//build tree from char freqency list
		H_tree * tree;
		buildTree(tree,charCount);

		//print header
		// printCharCount(charCount,ofs);

		map<unsigned char, string> codeMap;
		mapCodes(tree,codeMap,"");

		// printEncoding(ofs,codeMap);

		// ifs.seekg(0);

		// encodeChars(ifs,ofs,codeMap);

		bool complete = false;
		queue<size_t> bitStream;

		// printMap(codeMap);

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

bool decompression(string inFile, string outFile){
	fstream ifs(inFile, ios::in | ios::binary | ios::ate);
	fstream ofs(outFile, ios::out | ios::binary | ios::ate);

	if(ifs.is_open() && ofs.is_open()){
		ifs.seekg(ios_base::beg);
		ofs.seekp(ios_base::beg);

		char buff[1];
		tail_type extraBits;
		ifs.read( buff, sizeof(tail_type) );
		extraBits = (tail_type)(*buff);

		//get char frequency from compressed file
		ifs.seekg( ios_base::beg + sizeof(tail_type) );

		map<unsigned char, string> codeMap;
		readEncoding(ifs,codeMap);
		// string path;
		// mapCodes(tree,codeMap,"");
		// cout << "mapped codes" << endl;

		//reverse the map
		map<string,unsigned char> charMap;
		for(auto iter = codeMap.begin(); iter != codeMap.end(); iter++){
			charMap[iter->second] = iter->first;
		}
		
		// cout << "mapped chars" << endl;

		// printMap(codeMap);
		// printMap(charMap);

		bool complete = false;
		queue<unsigned char> byteStream;

		// cout << "done setting up for decoding" << endl;
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

bool mapFile(string inFile){
	fstream ifs(inFile, ios::in | ios::binary | ios::ate);

	if(ifs.is_open()){
		ifs.seekg(ios_base::beg);

		//get char frequency from compressed file
		freq_type charCount[CHARS] = {};
		countChars(charCount,ifs);

		//build tree from char frequency list
		H_tree * tree;
		buildTree(tree,charCount);

		map<unsigned char, string> codeMap;
		mapCodes(tree,codeMap,"");

		printMap(codeMap);
		return true;
	}
	else {
		return false;
	}
	return true;
}

void inputError(){
	cout << "Incorrect usage.\n";
	cout << "Run the program using the following options.\n";
	cout << "  Compression:\n";
	cout << "   -c SOURCE DESTINATION\n";
	cout << "  Decompression:\n";
	cout << "   -d SOURCE DESTINATION\n";
	cout << "  Print Map:\n";
	cout << "   -m SOURCE\n";
	cout << endl;
}

int main(int argc, char** argv){
	// thread test(printStuff,"hello world");
	// test.join();
	// // cout << "argc: " << argc << endl;
	// return 0;

	if(argc == 4){
		string option(argv[1]);
		string ifname(argv[2]);
		string ofname(argv[3]);

		if(option == "-c"){
			if(compression(ifname,ofname)) cout << "Done compressing file." << endl;
			else cout << "Unknown error when compressing." << endl;
		} 
		else if(option == "-d"){
			if(decompression(ifname,ofname)) cout << "Done decompressing file." << endl;
			else cout << "Unknown error when decomressing." << endl;
		} 
		else {
			inputError();
		}
	} else if(argc == 3){
		string option(argv[1]);
		string ifname(argv[2]);

		if(option == "-m"){
			if(mapFile(ifname)) cout << "Done mapping file." << endl;
			else cout << "Unknown error when mapping file." << endl;
		} 
		else {
			inputError();
		}
	} 
	else {
		inputError();
	}


		

 	// freq_type charCount[CHARS] = {};

	// string ifname,ofname;
	// fstream fs;

	// cout << "Enter an input filename: ";
	// cin >> ifname;
	// cout << "Enter an output filename: ";
	// cin >> ofname;

	

	
	// fs.open(fname, ios::in | ios::binary);

	// if(fs.is_open()){
	// 	cout << "opened file" << endl;
	// 	cout << "counting chars..";
	// 	countChars(charCount,fs);
	// 	cout << "...done" << endl;
	// 	cout << "creating tree..";
	// 	H_tree * tree;
	// 	buildTree(tree,charCount);
	// 	cout << "...done" << endl;
	// 	// cout << "printing inorder trasversal" << endl;
	// 	// vector<unsigned char> chars;
	// 	// inorder(tree,chars);
	// 	// for(auto i : chars){
	// 	// 	cout << (int)chars[i] << endl;
	// 	// }

	// 	map<unsigned char, string> codeMap;
	// 	string path;
	// 	cout << "mapping codes..";
	// 	cout.flush();
	// 	mapCodes(tree,codeMap,path);
	// 	cout << "...done" << endl;
	// 	for(auto i : codeMap){
	// 		cout << (size_t)i.first << " " << i.second << endl;
	// 	}



	// } else {
	// 	cout << "unable to open file" << endl;
	// }
	
	// int countChars = 0;
	// int countBits = 0;
	// for(int i = 0; i < CHARS; i++){
	// 	countBits += charCount[i];
	// 	// cout << "[" << charCount[i] << "] ";
	// 	if(charCount[i] != 0) countChars++;
	// }
	// cout << endl;
	// cout << "Total bits: " << countBits << endl;
	// cout << "Total chars: " << countChars << endl;


	return 0;
}