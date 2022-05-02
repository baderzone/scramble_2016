#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <time.h>

using namespace std;
void usage ()
{
        cout<<"Usage: duplication <read1.fq> <read2.fq> <clean1.fq> <clean2.fq> <dup.stat>\n";
		cout<<"Note:read1 length + read2 length < 225\n";
        exit(1);
}

typedef unsigned long long bit64_t; // 64-bit
typedef unsigned int bit8_t; //8-bit
typedef struct read1
{
	bit64_t AA[7];
}Read1;

struct Compare
{
	bool operator()(const read1 &a, const read1 &b)
	{
		for (int i=0;i<7 ;i++ )
		{
			if (a.AA[i] < b.AA[i])
			{
				return true;
			}
			else if (a.AA[i] > b.AA[i])
			{
				return false;
			}
		}
		return false;
	}
};


Read1 in_read1(string &seq,int n)
{
	Read1 base;
	base.AA[0]=0;
	base.AA[1]=0;
	base.AA[2]=0;
	base.AA[3]=0;
	base.AA[4]=0;
	base.AA[5]=0;
	base.AA[6]=0;
	int x=0;
	for (int m=0;m<n ;m++ )
	{
		x=m+1;
		for (int i=32*m;i<32*x ;i++ )
		{
			int k=((seq[i]&0x06)>>1);
			base.AA[m]=base.AA[m]<<2;
			base.AA[m]+=k;
		}
	}
	for (int i=32*n;i<seq.size() ;i++ )
	{
		int k=((seq[i]&0x06)>>1);
		base.AA[n]=base.AA[n]<<2;
		base.AA[n]+=k;
	}
	return base;
}

int main (int argc, char *argv[])
{
	if (argc<5) usage();
	string in_seq1=argv[optind++];
	string in_seq2=argv[optind++];
	string out_seq1=argv[optind++];
	string out_seq2=argv[optind++];
	string dup_stat=argv[optind++];

	string id1="", seq1="", qid1="", q1="", id2="", seq2="", qid2="", q2="";
	int read1_len=0, read2_len=0, total_len=0;

	ifstream infile1( in_seq1.c_str());
	if ( ! infile1)
	{
		cerr<<"fail to open input file"<<in_seq1<<endl;
	}
	getline ( infile1, id1, '\n');
	getline ( infile1, seq1, '\n');
	read1_len=seq1.length();

	ifstream infile2( in_seq2.c_str());
	if ( ! infile2)
	{
		cerr<<"fail to open input file"<<in_seq2<<endl;
	}
	getline ( infile2, id2, '\n');
	getline ( infile2, seq2, '\n');
	read2_len=seq2.length();

	total_len=read1_len+read2_len;
	if (total_len>224)
	{
		cerr<<"max pair reads length is 224"<<endl;
		exit(1);
	}
	if (total_len<1)
	{
		cerr<<"format error"<<endl;
		exit(1);
	}
	float size_len=total_len/32.0;
	cout<<"size_len"<<size_len<<endl;
	int size_num=(int)size_len;
	if (size_num==size_len)
	{
		size_num--;
	}
	cout<<"size_num:"<<size_num<<endl;
	infile1.close();
	infile2.close();
	id1="", id2="", seq1="", seq2="";

	infile1.open ( in_seq1.c_str());
	infile2.open ( in_seq2.c_str());

	ofstream outfile1 (out_seq1.c_str());
	if (! outfile1)
	{
		cerr << "fail to create output file" <<out_seq1<< endl;
	}
	ofstream outfile2 (out_seq2.c_str());
	if (! outfile2)
	{
		cerr << "fail to create output file" <<out_seq2<< endl;
	}
	ofstream stat (dup_stat.c_str());
	if (! stat)
	{
		cerr << "fail to create stat file" <<dup_stat<< endl;
	}
	
	map<Read1, int, Compare> read_count;
	int clean=0,total_num=0,dup=0;
	time_t start_time,end_time;
	time(&start_time);
	while ( getline ( infile1, id1, '\n') )
	{
		if (id1[0] == '@')
		{
			getline ( infile1, seq1, '\n');
			getline ( infile1, qid1, '\n');
			getline ( infile1, q1, '\n');
			getline ( infile2, id2, '\n');
			getline ( infile2, seq2, '\n');
			getline ( infile2, qid2, '\n');
			getline ( infile2, q2, '\n');
			total_num++;
			string read=seq1+seq2;
			Read1 read_seq=in_read1(read,size_num);
			//cout<<read_count[read_seq]<<endl;
			if (read_count[read_seq]>0)
			{
				read_count[read_seq]++;
				//cout<<id1<<endl<<seq1<<endl<<qid1<<endl<<q1<<endl;
				//cout<<id2<<endl<<seq2<<endl<<qid2<<endl<<q2<<endl;
			}
			else
			{
				read_count[read_seq]=1;
				clean++;
				outfile1<<id1<<endl<<seq1<<endl<<qid1<<endl<<q1<<endl;
				outfile2<<id2<<endl<<seq2<<endl<<qid2<<endl<<q2<<endl;
			}
		}
	}
	map<Read1, int>::const_iterator map_it = read_count.begin();
	while (map_it != read_count.end())
	{
		if (map_it->second >1)
		{
			dup+=map_it->second;
		}
		//cout<<map_it->second<<endl;
		map_it++;
	}
	time(&end_time);
	stat<<"Total_reads:"<<total_num<<endl<<"Duplicate_reads:"<<dup<<endl<<"Clean_reads:"<<clean<<endl<<"Used_time:"<<end_time-start_time<<endl;
	
}

