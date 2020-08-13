#include <iostream>
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
#include <thread>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
using namespace std;
#ifdef WIN32
#include <direct.h>;
#else
#include <sys/stat.h>;
#endif

void createDir(const char* dir){
	int r;
#ifdef WIN32

	_mkdir(dir);
#else
	r = mkdir(dir, 0755);
#endif

}

static void show_usage(std::string name)
{
	std::cerr << "VCF to BEDp Converter. v.0.1\n";

	std::cerr << "Usage: " << name << " <option(s)>\n"

			<< "Options:\n"
			<< "\t-h,--help                               Show this help message\n"
			<< "\t-i,--vcf <INPUT FILE>                   Input file path in compressed VCF format (.gz)\n"
			<< "\t-o,--output_file <OUTPUT PREFIX>        Output prefix\n"
			<< "\t-o,--output_folder <OUTPUT DIRECTORY>   Output folder\n"

			<< std::endl;
}


int main(int argc, char* argv[]){

	char* output_folder;
	char* input_vcf;
	char* output_prefix;

	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "--help")) {
			show_usage(argv[0]);
			return 0;
		} else if ((arg == "-i") || (arg == "--input")) {
			if (i + 1 < argc) {
				input_vcf = argv[++i];
			}
		}
		else if ((arg == "-O") || (arg == "--output_folder")) {
			if (i + 1 < argc) {
				output_folder = argv[++i];
			}
		}

		else if ((arg == "-o") || (arg == "--output_file")) {
			if (i + 1 < argc) {
				output_prefix = argv[++i];
			}
		}


	}

	if (argc < 3){
		show_usage(argv[0]);

		return -1;
	}
	createDir(output_folder);

	clock_t start = clock();
	std::ifstream input_vcf_file(input_vcf, std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
	inbuf.push(boost::iostreams::gzip_decompressor());
	inbuf.push(input_vcf_file);
	std::istream vcf(&inbuf);
	stringstream ss1;
	ss1 << output_folder << "/" << output_prefix;
	ss1 << ".fam";

	stringstream ss2;
	ss2 << output_folder << "/" << output_prefix;
	ss2 << ".bim";

	stringstream ss3;
	ss3 << output_folder << "/"  << output_prefix;
	ss3 << ".bed";

	std::ofstream fam(ss1.str());
	std::ofstream bim(ss2.str());
	std::ofstream bedp(ss3.str(), ios::binary|ios::out);


	std::bitset<8> x;
	x.reset();
	x.set(2);
	x.set(3);
	x.set(5);
	x.set(6);
	bedp<<static_cast<uint_fast8_t>(x.to_ulong());
	x.reset();
	x.set(0);
	x.set(1);
	x.set(3);
	x.set(4);
	bedp<<static_cast<uint_fast8_t>(x.to_ulong());
	x.reset();
	x.set(0);
	bedp<<static_cast<uint_fast8_t>(x.to_ulong());


	const char delim = '\t';
	int header = 0;
	int lineCount = 0;
	float leftOver = 0;
	std::string line;
	while(std::getline(vcf, line)){
		std::stringstream s(line);
		if (line.substr(0,2) == "##"){
			header = header + 1;
		}
		else{
			if (line.substr(0,1) == "#"){
				int count = 0;
				std::string token;
				while(std::getline(s,token,delim)){
					if (count > 8){
						fam<<token<<"\t"<<token<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"-9"<<"\n";
					}
					count = count + 1;
				}
				fam.close();
				leftOver = float(count - 9) / 4.0;
				leftOver = leftOver - int(leftOver);
			}
			else{
				std::string token2;
				vector<string> site;
				while(std::getline(s,token2,delim)){
					site.push_back(token2);
				}
				if (site[8] != "GT"){
					cout<<site[1]<<"\t"<<lineCount<<"\n";
					continue;
				}
				bim<<site[0]<<"\t"<<site[2]<<"\t"<<"0"<<"\t"<<site[1]<<"\t"<<site[4]<<"\t"<<site[3]<<"\n";
				for (int x = 9; x < site.size(); x = x + 4){
					std::bitset<8> a;
					std::bitset<8> b;
					std::bitset<8> c;
					std::bitset<8> d;
					a.reset();
					b.reset();
					c.reset();
					d.reset();
					if (site[x] == "0|0"){
						a.reset();
					}
					else if (site[x] == "0|1"){
						a.set(0);
					}
					else if (site[x] == "1|1"){
						a.set(0);
						a.set(1);
					}
					else if (site[x] == "1|0"){
						a.set(1);
					}
					if (site[x + 1] == "0|0"){
						b.reset();
					}
					else if (site[x + 1] == "0|1"){
						b.set(2);
					}
					else if (site[x + 1] == "1|1"){
						b.set(2);
						b.set(3);
					}
					else if (site[x + 1] == "1|0"){
						b.set(3);
					}
					if (site[x + 2] == "0|0"){
						c.reset();
					}
					else if (site[x + 2] == "0|1"){
						c.set(4);
					}
					else if (site[x + 2] == "1|1"){
						c.set(4);
						c.set(5);
					}
					else if (site[x + 2] == "1|0"){
						c.set(5);
					}
					if (site[x + 3] == "0|0"){
						d.reset();
					}
					else if (site[x + 3] == "0|1"){
						d.set(6);
					}
					else if (site[x + 3] == "1|1"){
						d.set(6);
						d.set(7);
					}
					else if (site[x + 3] == "1|0"){
						d.set(7);
					}
					std::bitset<8> byte = (a | b | c | d);
					bedp<<static_cast<uint_fast8_t>(byte.to_ulong());
				}
				if (leftOver == 0.25){
					std::bitset<8> a;
					std::bitset<8> b;
					std::bitset<8> c;
					std::bitset<8> d;
					a.reset();
					b.reset();
					c.reset();
					d.reset();
					if (site[site.size()-1] == "0|0"){
						a.reset();
					}
					else if (site[site.size()-1] == "1|1"){
						a.set(0);
						a.set(1);
					}
					else if (site[site.size()-1] == "0|1"){
						a.set(0);
					}
					else if (site[site.size()-1] == "1|0"){
						a.set(1);
					}
					std::bitset<8> byte2 = (a | b| c | d);
					bedp<<static_cast<uint_fast8_t>(byte2.to_ulong());
				}
				if (leftOver == .5){
					std::bitset<8> a;
					std::bitset<8> b;
					std::bitset<8> c;
					std::bitset<8> d;
					a.reset();
					b.reset();
					c.reset();
					d.reset();
					if (site[site.size()-2] == "0|0"){
						a.reset();
					}
					else if (site[site.size()-2] == "1|1"){
						a.set(0);
						a.set(1);
					}
					else if (site[site.size()-2] == "0|1"){
						a.set(0);
					}
					else if (site[site.size()-2] == "1|0"){
						a.set(1);
					}
					if (site[site.size()-1] == "0|0"){
						b.reset();
					}
					else if (site[site.size()-1] == "1|1"){
						b.set(2);
						b.set(3);
					}
					else if (site[site.size()-1] == "0|1"){
						b.set(2);
					}
					else if (site[site.size()-1] == "1|0"){
						b.set(3);
					}
					std::bitset<8> byte2 = (a | b| c | d);
					bedp<<static_cast<uint_fast8_t>(byte2.to_ulong());
				}
				if (leftOver == .75){
					std::bitset<8> a;
					std::bitset<8> b;
					std::bitset<8> c;
					std::bitset<8> d;
					a.reset();
					b.reset();
					c.reset();
					d.reset();
					if (site[site.size()-3] == "0|0"){
						a.reset();
					}
					else if (site[site.size()-3] == "1|1"){
						a.set(0);
						a.set(1);
					}
					else if (site[site.size()-3] == "0|1"){
						a.set(0);
					}
					else if (site[site.size()-3] == "1|0"){
						a.set(1);
					}
					if (site[site.size()-2] == "0|0"){
						b.reset();
					}
					else if (site[site.size()-2] == "1|1"){
						b.set(2);
						b.set(3);
					}
					else if (site[site.size()-2] == "0|1"){
						b.set(2);
					}
					else if (site[site.size()-2] == "1|0"){
						b.set(3);
					}
					if (site[site.size()-1] == "0|0"){
						c.reset();
					}
					else if (site[site.size()-1] == "1|1"){
						c.set(4);
						c.set(5);
					}
					else if (site[site.size()-1] == "0|1"){
						c.set(4);
					}
					else if (site[site.size()-1] == "1|0"){
						c.set(5);
					}
					std::bitset<8> byte2 = (a | b| c | d);
					bedp<<static_cast<uint_fast8_t>(byte2.to_ulong());
				}
			}
			lineCount = lineCount + 1;
		}
	}
	input_vcf_file.close();
	bim.close();
	bedp.close();
	clock_t end = clock();
	double time = ((double) (end-start)) / CLOCKS_PER_SEC;
	cout<<time<<endl;
}
