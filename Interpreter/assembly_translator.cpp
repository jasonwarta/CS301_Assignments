#include <cstdio>
#include <cstdlib>
#include <ctype.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "lib/inc.h"

using std::ios;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;

//initialzie global variables for text and filesize
unsigned char * src;
size_t size = 0;

// Allocate writeable, executable memory:
#include <sys/mman.h>
unsigned char *bytes=(unsigned char *)mmap(0,16384,
		PROT_READ+PROT_WRITE+PROT_EXEC,
		MAP_ANON+MAP_SHARED,
		-1,0);

int foo(int input) {
	//move text into new array to be called
	for (unsigned int i=0;i<size;i++) bytes[i]=src[i];

	// Run the new function
	typedef long  (*fnptr)  (int);
	fnptr f=(fnptr)bytes;
	return f(input);
}

int main(){
	//prompt for and store filename
	cout << "Enter filename: ";
	string fname;
	cin >> fname;

	//prompt for int input
	char toggle;
	int input;
	cout << "Do you want to pass an int to the program? (y/n): ";
	cin >> toggle;
	//store input if the user provides it
	if(toggle == 'y' || toggle == 'Y'){
		cout << "Enter the int: ";
		cin >> input;
	}

	//open file
	ifstream ifs(fname, ios::binary);

	unsigned char * text = 0;
	size;

	if(ifs.is_open()){
		//get filesize for intizalizing array
		ifs.seekg(0,ios::end);
		size = ifs.tellg();

		cout << "file size: " << size << " bytes" << endl;
		ifs.seekg(0,ios::beg);

		//move contents of file into array of characters
		text = new unsigned char[size];
		ifs.read( (char*)text,size );

		ifs.close();

		cout << "finished reading file" << endl;

		src = new unsigned char[size];

		for(int i = 0; i < size; i++){
			src[i] = text[i];
		}

		cout << "done moving code to array. calling code." << endl;

		//call loaded program and print output
		cout << foo(input) << endl;

	} else {//file didn't open correctly
		cout << "unable to open file." << endl;
	}

	return 0;
}

