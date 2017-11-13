#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <array>
#include <queue>
#include <vector>
#include <map>
#include <list>
#include <stack>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include "chartobits.h"

// cd /mnt/c/Users/AJ\ Steenhoek/Documents/Fall_2017/Data_Structures/lab10
// g++ -std=c++14 huffman.cpp -ohuffman
using namespace std;

namespace Huffman {
#define ZIPEXT ".zip"
#define TREEEXT ".tree"
	
struct Info { unsigned char c; int freq; };

class Node {
public:
	unsigned char c;
	unsigned int freq;
	string code;
	Node *left, *right;
	Node(): freq(-1), left(nullptr), right(nullptr), code("") { ; }
	bool isLeaf() { return left == nullptr && right == nullptr; }
}; // end of struct Node

void inorderPrint(Node *root) {
    if (root)
    {
        inorderPrint(root->left);
        printf("%c ", root->c);
        inorderPrint(root->right);
    }
}

string inorderEncode(Node *root) {
	static string codestring;
	if(!root->isLeaf()) {
		codestring.append(1, '0');
		inorderEncode(root->left);
		inorderEncode(root->right);
	}
	else {
		codestring.append(1, '1');
		codestring.append(1, root->c);
	}
	return codestring;
}

struct CompressionTree {
	Node* root;
}; // end of struct CompressionTree

class Greater {
public:
	bool operator() (Node* node1, Node* node2) {
		return node1->freq > node2->freq;
	}
}; // end of struct Greater

std::priority_queue<Node*, std::vector<Node*>, Greater> minheap;
std::map<char, std::string> codes;
CompressionTree tree;
//string bits;

// An iterative process to encode all the chars in the file
void encode(Node* node) {
	static std::string bits = "";
	if (node->right != NULL)
	{
		bits += "1";
		// std::cout << "parent: " << node->freq << std::endl;
		encode(node->right);
		bits = bits.substr(0, bits.size() - 1);
	}
	if (node->left != NULL)
	{
		bits += "0";
		// std::cout << "parent: " << node->freq << std::endl;
		encode(node->left);
		bits = bits.substr(0, bits.size() - 1);
	}
	if(!node->left && !node->right)
	{
		codes[node->c] = bits;
		// std::cout << node->c << " \n";
	}
}

void printCodes() {
	for(auto i : codes) {
		cout << i.first << ": " << i.second << endl;
	}
}

string bitString(const char* filename) {

	std::ifstream ifs (filename, std::ifstream::in);
	
	char c = ifs.get();
	string bits = "";

	while (ifs.good()) {
	bits += codes[c];
	//std::cout << c;
	c = ifs.get();
	}
	bits += codes[c];

	//cout << bits.size() << endl;
	//tag on the padding bits
	int rem = 8 - bits.size() % 8;
	//cout << rem << endl;
	std::bitset<8> header(rem);

	//cout << 8 + rem << endl;

	string pad = "";
	while(rem--) {
		pad += "1";
	}
	
	bits = header.to_string() + pad + bits;

	ifs.close();
	//cout << bits.size() << endl;
	return bits;
}

void printFile(const char* filename) {
	ifstream file(filename);
	string contents;
	file >> contents;
	cout << contents << endl;
}

void buildTree(std::array<int, 256> &counts) {
	//cout << "buildTree" << endl;
	while(!minheap.empty())
		minheap.pop();
	tree.root = nullptr;
	for (int c=0; c < 256; ++c) {
		if (counts[c] > 0) {
			Node* node = new Node();
			node->c = c;
			node->freq = counts[c];
			minheap.push(node);
		}
	}

	while (minheap.size() > 1) {
		Node *node = new Node();
		Node *leftchild = minheap.top();
		minheap.pop();
		Node *rightchild = minheap.top();
		minheap.pop();

		node->left = leftchild;
		node->right = rightchild;
		node->freq = leftchild->freq + rightchild->freq;
		node->c = -1;

		minheap.push(node);
	}

	tree.root = minheap.top();
}

// read chars from file into separate counts array
bool readCounts(const char *filename, std::array<int, 256> &counts) {
	FILE *file = fopen(filename, "rb");
	if (nullptr == file) { return false; }

	for (int &cnt : counts) { cnt = 0; }

	Info info;
	while (1 == fread(&info, sizeof(Info), 1, file)) {
		counts[info.c] = info.freq;
	}

	fclose(file);
	return true;
}

// NOTE: pops everything off of the heap
void printHeap() {
	while (!minheap.empty()) {
		printf("%c ", minheap.top()->c);
		minheap.pop();
	}
	printf("\n");
}

bool writeCounts(const char *filename, const std::array<int, 256> &counts) {
	FILE *file = fopen(filename, "wb");
	if (nullptr == file) { return false; }

	Info info;
	for (int c=0; c < 256; ++c) {
		if (counts[c] > 0) {
			info.c = (unsigned char)c;
			info.freq = counts[c];
			fwrite(&info, sizeof(Info), 1, file);
		}
	}

	fclose(file);
	
	return true;
}

bool getCounts(const char *filename, std::array<int, 256> &counts) {
	FILE *file = fopen(filename, "rb");
	if (nullptr == file) { return false; }
	
	for (int &cnt : counts) { cnt = 0; }

	unsigned char c;
	while (1 == fread(&c, 1, 1, file)) {
		++counts[c];
	}

	fclose(file);

	return true;
}

//fail
bool readBits(string filename, string &contents) {
	//inorderPrint(tree.root);
	//printCodes();
	FILE *file = fopen(filename.c_str(), "rb");
	if (nullptr == file) { return false; }
	unsigned char skipbits;
	fread(&skipbits, sizeof(unsigned char), 1, file);

	queue<bool> bits;
	unsigned char c;
	unsigned char shifter = 128;
	while (1 == fread(&c, sizeof(unsigned char), 1, file)) {
		while (shifter > 0) {
		bool bit = c & shifter;
		shifter >>= 1;
		bits.push(bit);
		}
		shifter = 128;
	}
	bits.push(0);
	string outfilename = filename.substr(0, filename.size() - 4);
	while (skipbits--) {
		bits.pop();
	}

	// while(!bits.empty()) {
	// 	cout << bits.front();
	// 	bits.pop();
	// }
	//cout << endl;
	Node *node = tree.root;
	//cout << "decode file" << endl;
	while (!bits.empty()) {
		if(node->isLeaf()) {
			//FILE *outfile = fopen(outfilename.c_str(), "wb");
			//fwrite(&(node->c), 1, 1, outfile);
			//fclose(outfile);
			contents += string(1,node->c);
			// cout << node->c << endl;
			// add char to file
			node = tree.root;
			continue;
		}
		else {
			if (bits.front() == 1) {
				node = node->right;
			}
			else {
				node = node->left;
			}
		}
		bits.pop();
	}
	ofstream outfile(outfilename);
	outfile << contents;
	outfile.close();
	
	fclose(file);

	return false;
}

} // end of namespace Huffman

//using namespace std;

int main(int argc, const char *argv[]) {
	if (argc < 3) { 
		puts("usage/compress: huffman -c fileToCompress\nusage/decompress: huffman -d fileToDecompress"); 
		exit(1); 
	}
	//cout << "main" << endl;
	if (strcmp(argv[1], "-d") == 0) {
		
		string unzipFilename(argv[2]);
		string treeFilename(argv[2]);
		for(int i = 0; i < 4; i++)
			treeFilename.pop_back();
		treeFilename += TREEEXT;

		std::array<int, 256> counts;
		
		Huffman::readCounts(treeFilename.c_str(), counts);
		// for (int k=0; k < 256; ++k) {
		// 	if (counts[k] > 0) {
		// 		printf("%3d %c: %d\r\n", k, ((32 < k && k < 127) ? ((char)k) : ' '), counts[k]);
		// 	}
		// }
		Huffman::buildTree(counts);
		string contents = "";
		Huffman::readBits(unzipFilename.c_str(), contents);
		
		//cout << contents << endl;
		// Huffman::inorderPrint(Huffman::tree.root);
		// cout << endl;
	} 
	else if (strcmp(argv[1], "-c") == 0) {
		std::array<int, 256> counts;
		if (!Huffman::getCounts(argv[2], counts)) {
			puts("The file could not be opened.");
			exit(2);
		}
		Huffman::buildTree(counts);
		// for (int k=0; k < 256; ++k) {
		// 	if (counts[k] > 0) {
		// 		printf("%3d %c: %d\r\n", k, ((32 < k && k < 127) ? ((char)k) : ' '), counts[k]);
		// 	}
		// }

		std::string treeFilename(argv[2]);
		std::string zipFilename(argv[2]);
		treeFilename += TREEEXT;
		zipFilename += ZIPEXT;
		// writes freq to .tree file
		Huffman::writeCounts(treeFilename.c_str(), counts);
		Huffman::encode(Huffman::tree.root);
		// Huffman::printCodes();
		// Huffman::inorderPrint(Huffman::tree.root);
		// cout << endl;
		string bitString = Huffman::bitString(argv[2]);
		// cout << bitString << endl;
		// Huffman::inorderPrint(Huffman::tree.root);
		// cout << endl;
		CharToBits ch(bitString);
		int len = ch.insertBits(zipFilename.c_str());
		// if(len == bitString.size()) {
		// 	cout << "Compression successful" << endl;
		// }
		// else
		//  cout << "Compression failed" << endl;
	} else {
		puts("usage/compress: huffman -c fileToCompress\nusage/decompress: huffman -d fileToDecompress"); 
		exit(1); 
	}
	return 0;
}
