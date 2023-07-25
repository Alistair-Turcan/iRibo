#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include <cstring>
#include <memory>
#include <random>
#include <algorithm>
#include <set>
#include <string_view>
#include <thread>
#include <chrono>
#include <omp.h>
#include <array>
#include <filesystem>

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

//constexpr string_view GENOME_ANNOTATION_PATH {"GCF_000002035.6_GRCz11_genomic.gtf"};
//constexpr string_view GENOME_PATH {"GCF_000002035.6_GRCz11_genomic.fna"};

//Four tasks of iRibo
//1. construct candidate ORFs from genome and transcript annotation
//2. read in ribo-seq mapping data and associate with genome 
//3. remap to p-site based on canonical genes
//4. map reads to ORFs
//5. assess candidate translation status

class GTF
{
public:
	string contig;
	string source;
	int start;
	int end;
	int strand;
	int exon_number = 0;
	int chr = -1;
	string annotation_type;
	string gene_id;
	string parent_id;
	string transcript_id;
	string gene_biotype;
	string ccds_id;
	vector<int> read_count;
	bool lifted = false;
	bool spliced_gene=false;
	bool operator<(const GTF& ann)
	{
		if (ann.chr == this->chr)
		{
			return this->start < ann.start;
		}
		return this->chr < ann.chr;
	}
	GTF()
	{
		read_count = vector<int>(3);
	}
};

class Exon
{
public:
	int start=-1;
	int end=-1;
	string seq;
	bool is_coding = false;
	int num_agree = 0;
	int pos_agree = 0;
	Exon()
	{
		//seq = vector<vector<int>>(2);
	}
	bool operator<(const Exon& exon)
	{
		return this->start < exon.start;
	}

};

class Transcript
{
public:
	int chr = -1;
	string transcript_id;
	int strand = -1;
	vector<Exon> exons;
	string seq;
	bool bad = false;
};

class ORF
{
public:
	int start=-1;
	int end=-1;
	int rfc =0;
	vector<int> frame;
};

class Read_Coord
{
public:
	int length1;
	int chr;
	int start;
};

class CandidateORF
{
public:
	int chr = -1;
	string chr_str;
	string transcript_id;
	string gene_id;
	string secondary_gene_id;
	int start_codon_pos=-1;
	int stop_codon_pos = -1;

	int strand = -1;
	vector<Exon> exons;
	int num_agree = 0;
	int pos_agree = 0;
	string species="none";
	bool intersect_coding = false;
	bool lncrna = false;
	bool protein_coding = false;
	bool processed_pseudogene = false;
	bool unprocessed_pseudogene = false;
	string biotype_other = "X";
	string intersect_gene = "X";
	string intersect_five_utr = "X";
	string intersect_three_utr = "X";
	string antisense_gene = "X";
	string uorf_id = "X";
	string dorf_id = "X";
	string CDS_intersect = "X";
	string CDS_biotype = "X";
	int orf_length = 0;
	int age_class = -1;
	string orf_annotation_type;
	string orf_gene_biotype;

	CandidateORF()
	{

	}
	bool operator<(const CandidateORF& gene)
	{
		
		if (gene.chr == this->chr)
		{
			
			return this->start_codon_pos < gene.start_codon_pos;
		}
		return this->chr < gene.chr;
	}
};

class GeneModel
{
public:
	int chr = -1;
	//string chr_str;
	//string transcript_id;
	string gene_id;
	//string secondary_gene_id;
	int start_codon_pos=-1;
	int stop_codon_pos = -1;
	//vector<int> exon_starts;
	//vector<int> exon_stops;
	int strand = -1;
	vector<int> frame;
	//vector<int> read_count;
	vector<int> pos_read_count;
	//vector<int> read_priority;
	//vector<int> read_priority_fake;
	//vector<int> read_priority_reads;
	//vector<int> read_priority_reads_fake;
	int first_max=0;
	int first_max_scrambled=0;
	int total_frames=0;
	int total_frames_scrambled=0;
	vector<int> frame_reads;
	vector<int> frame_reads_scrambled;

	vector<Exon> exons;
	//int num_agree = 0;
	//int pos_agree = 0;
	//string species="none";
	//bool intersect_coding = false;
	//bool lncrna = false;
	//bool protein_coding = false;
	//bool processed_pseudogene = false;
	//bool unprocessed_pseudogene = false;
	//string biotype_other = "X";
	//string intersect_gene = "X";
	//string intersect_five_utr = "X";
	//string intersect_three_utr = "X";
	//string antisense_gene = "X";
	//string uorf_id = "X";
	//string dorf_id = "X";
	//string CDS_intersect = "X";
	//string CDS_biotype = "X";
	//int orf_length = 0;
	//int age_class = -1;
	//string orf_annotation_type;
	//string orf_gene_biotype;
	//bool has_start_reads = false;
	//int first_half_reads = 0;
	//int second_half_reads=0;
	//int first_half_reads_fake = 0;
	//int second_half_reads_fake=0;	
	//int seq_len=0;
	double p_value;
	double p_value_scrambled;
	
	GeneModel()
	{
		//read_count = vector<int>(3);
		//read_priority = vector<int>(3);
		//read_priority_fake = vector<int>(3);
		//read_priority_reads = vector<int>(3);
		//read_priority_reads_fake = vector<int>(3);
		//max_codon = vector<int>(3);
		//max_codon_scrambled = vector<int>(3);
		frame_reads = vector<int>(3);
		frame_reads_scrambled = vector<int>(3);

		//seq = vector<vector<int>>(SPECIES_COUNT);
		//exon_seqs = vector<vector<vector<int>>>(2);
	}
	bool operator<(const GeneModel& gene)
	{
		if (gene.chr == this->chr)
		{
			return this->start_codon_pos < gene.start_codon_pos;
		}
		return this->chr < gene.chr;
	}
};

class CDS
{
public:
	int start=-1;
	int end=-1;
	map<int,int> exon_starts;
	map<int,int> exon_stops;
};

class Read
{
public:
	int chr;
	int start;
	int length;
	int strand = -1;
	int flag = -1;
	int count = 1;
};

class Study
{
public:
	string path;
	//string srp;
	//string srr;
	//vector<string> files;
	//bool preprocessed;
	//string linker;
	//vector<map<int, int>> frames;
	//vector<int> file_pass;
	//bool wildtype = true;
	//bool control = false;
};



constexpr std::array<char, 128> complement = [] {
    std::array<char, 128> lookup{};
    for (auto & c : lookup) { c = 'N'; }
    lookup['A'] = 'T';
    lookup['C'] = 'G';
    lookup['G'] = 'C';
    lookup['T'] = 'A';
    return lookup;
}();

void reverse_complement(string& sequence) {
    string result(sequence);
    transform(sequence.rbegin(), sequence.rend(), result.begin(), [](char c) { return complement[c]; });
    sequence = result;
}

void reverse_complement_old(string &seq)
{
std::map<char, char> comp = {
    {'A', 'T'},
    {'C', 'G'},
    {'G', 'C'},
    {'T', 'A'},
    {'B', 'V'},
    {'D', 'H'},
    {'H', 'D'},
    {'K', 'M'},
    {'M', 'K'},
    {'N', 'N'},
    {'R', 'Y'},
    {'S', 'S'},
    {'V', 'B'},
    {'W', 'W'},
    {'Y', 'R'}
};
	string rc_seq;
	for (int i = 0; i < seq.size(); i++)
	{
		rc_seq += (comp[seq[i]]);
	}
	reverse(rc_seq.begin(), rc_seq.end());
	seq = rc_seq;
}
void reverse_complement(vector<int> &rc, const vector<int> &seq)
{
	vector<int> comp = { 3,2,1,0,4,5,6 };
	for (int i = seq.size() - 1; i >= 0; i--)
	{
		if (seq[i] > comp.size())
		{
			cout << "\nmismatch:" << seq.size() << " " << seq[i];
			//getchar();
		}
		else
		{
			rc.push_back(comp.at(seq[i]));
		}
	}
}

void reverse_complement(vector<int> &seq)
{
	vector<int> rc;
	reverse_complement(rc, seq);
	seq = rc;
}

bool is_integer(const string & s)
{
	if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

	char * p;
	int n = strtol(s.c_str(), &p, 10);

	return (*p == 0);
}

void split(const string &s, char delim, vector<string> &elems)
{
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) 
	{
		elems.push_back(item);
	}
}
unordered_map<string, string> parseArguments(int argc, char* argv[]) {
    unordered_map<string, string> arguments;

    for (int i = 1; i < argc; ++i) {
        string arg(argv[i]);

        if (arg.substr(0, 2) == "--") {
            string::size_type delimiterPos = arg.find('=');
            if (delimiterPos != string::npos) {
                string key = arg.substr(2, delimiterPos - 2);
                string value = arg.substr(delimiterPos + 1);
                arguments[key] = value;
            } else {
                string key = arg.substr(2);
                arguments[key] = "";
            }
        }
    }

    return arguments;
}
// Function to retrieve argument value with a default or throw an error if mandatory and not specified
string getArg(const std::unordered_map<std::string, std::string>& arguments,  const std::string& arg, const string& defaultValue, bool mandatory) {
    auto it = arguments.find(arg);
    if (it != arguments.end()) {
        // Argument found, set the output variable to its value
		return it->second;
    } else {
        if (mandatory) {
            // Argument not found and mandatory, throw an error
            throw std::runtime_error("Argument " + arg + " is required but not specified.");
			exit(0);
        } else {
            // Argument not found, set the output variable to the default value
			return defaultValue;
        }
    }
}

/*
void read_gff(vector<GTF>& gtfs, string filename, vector<string, int> &chr_labels)
{

	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		if (line[0] != '#')
		{
			vector<string> columns;
			split(line, '\t', columns);
			gtfs.push_back(GTF());
			if (!chr_labels.count(columns[0]))
			{
				continue;
			}
			gtfs.back().chr = chr_labels.at(columns[0]) - 1;
			gtfs.back().annotation_type = columns[2];
			if (columns[2] == "pseudogene")
			{
				gtfs.back().orf_classification = "pseudogene";
			}
			else if (columns[2] == "transposable_element_gene")
			{
				gtfs.back().orf_classification = "transposable_element_gene";
			}
			else if (columns[2] == "LTR_retrotransposon")
			{
				gtfs.back().orf_classification = "LTR_retrotransposon";
			}
			else if (columns[2] == "blocked_reading_frame")
			{
				gtfs.back().orf_classification = "blocked_reading_frame";
			}
			else if (columns[2] == "tRNA_gene")
			{
				gtfs.back().orf_classification = "tRNA_gene";
			}
			else if (columns[2] == "long_terminal_repeat")
			{
				gtfs.back().orf_classification = "long_terminal_repeat";
			}
			else if (columns[2] == "snoRNA_gene")
			{
				gtfs.back().orf_classification = "snoRNA_gene";
			}
			else if (columns[2] == "snRNA_gene")
			{
				gtfs.back().orf_classification = "snRNA_gene";
			}
			else if (columns[2] == "ARS")
			{
				gtfs.back().orf_classification = "ARS";
			}
			else if (columns[2] == "ncRNA_gene")
			{
				gtfs.back().orf_classification = "ncRNA_gene";
			}

			gtfs.back().start = stoi(columns[3]) - 1;
			gtfs.back().end = stoi(columns[4]) - 1;
			
			if (columns[6] == "+")
			{
				gtfs.back().strand = 0;
			}
			else
			{
				gtfs.back().strand = 1;
			}
			vector<string> info;
			split(columns[8], ';', info);
			for (int i = 0; i < info.size(); i++)
			{
				vector<string> components;
				split(info[i], '=', components);
				if (components[0] == "ID")
				{
					gtfs.back().ID = components[1];
				}
				if (components[0] == "Name")
				{
					gtfs.back().Name = components[1];
				}
				if (components[0] == "Parent")
				{
					gtfs.back().Parent = components[1];
				}
				if (components[0] == "orf_classification")
				{
					gtfs.back().orf_classification = components[1];
				}
			}
		}
	}
}
*/
string removeSuffix(const string& str, const string& suffix) {
    if (str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0) {
        return str.substr(0, str.size() - suffix.size());
    }
    return str;
}

void read_gff3(vector<GTF>& gtfs, string filename, map<string, int>& chr_labels)
{

	ifstream file(filename);
	string line;
	int gene_count = 0;
	bool prev_teg = false; //if previous element was transposable_element_gene
	bool prev_CDS = false;
	bool prev_Intron = false;
	while (getline(file, line))
	{
		if (line[0] != '#')
		{

			vector<string> columns;
			split(line, '\t', columns);
			if (!chr_labels.count(columns[0]))
			{
				continue;
			}
			
			//Handle introns by extending the start/stop of the gene
			if(columns[2]=="intron" && prev_CDS){
				prev_Intron=true;
				continue;
			}
			prev_CDS=false;
			
			if(columns[2]=="CDS" && prev_Intron){
				if(columns[6] == "+"){
					gtfs.back().end = stoi(columns[4]) - 1;
				} else {
					gtfs.back().start = stoi(columns[3]) - 1;
				}
				gtfs.back().spliced_gene=true;
				prev_Intron=false;
				prev_CDS=true;
				continue;
			}
			gtfs.push_back(GTF());
			gtfs.back().chr = chr_labels[columns[0]];
			gtfs.back().contig = columns[0];
			
			gtfs.back().source = columns[1];
			gtfs.back().annotation_type = columns[2];



			
			gtfs.back().start = stoi(columns[3]) - 1;
			gtfs.back().end = stoi(columns[4]) - 1;
			
			if (columns[6] == "+")
			{
				gtfs.back().strand = 0;
			}
			else
			{
				gtfs.back().strand = 1;
			}
			vector<string> info;
			split(columns[8], ';', info);
			string gene_id;
			bool is_gene = false;
			for (int i = 0; i < info.size(); i++)
			{
				vector<string> components;
				split(info[i], '=', components);
				if (components[0] == "ID")
				{
					//gtfs.back().gene_id = components[1];
					//gene_id = components[1];
				}
				if (components[0] == "gene")
				{
					//gtfs.back().gene_id = components[1];
					//is_gene = true;
					//gene_id = components[1];

				}
				if (components[0] == "Name")
				{
					//gtfs.back().Name = components[1];
					gene_id = components[1];
					if(prev_teg){
						gene_id = removeSuffix(gene_id, "_CDS");
						gtfs.back().gene_id = gene_id;
						gene_count++;
						prev_teg = false;
					};

				}
				if (components[0] == "Parent")
				{
					gtfs.back().transcript_id = components[1];
				}
				if (components[0] == "orf_classification" && gtfs.back().annotation_type == "CDS")
				{
					prev_CDS=true;
					if(components[1] == "Verified" || components[1] == "Uncharacterized"){
						gene_id = removeSuffix(gene_id, "_CDS");
						gtfs.back().gene_id = gene_id;
						//cout<<"\n" << gtfs.back().gene_id;
						gene_count++;
						//cout <<"\nFound a potential gene: " << to_string(gene_count);
					} else { //Deal with dubious, pseudogene, etc.
						gtfs.back().annotation_type = "";
					}
					
				}
				
			}
			if(gtfs.back().annotation_type == "transposable_element_gene"){
				prev_teg = true;
			}

		}
	}
}
void read_gtf_original(vector<GTF> &gtfs, string filename, map<string,int> &chr_labels, bool includes_chr_prefix)
{
	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		if (line[0] != '#')
		{
			vector<string> columns;
			split(line, '\t', columns);
			if(!includes_chr_prefix)
			{
				columns[0] = "chr"+columns[0];
			}
			if(!chr_labels.count(columns[0]))
			{
				continue;
			}
			gtfs.push_back(GTF());
			gtfs.back().chr = chr_labels.at(columns[0]);
			gtfs.back().contig = columns[0];
			
			gtfs.back().source = columns[1];
			gtfs.back().annotation_type = columns[2];
			gtfs.back().start = stoi(columns[3])-1;

			gtfs.back().end = stoi(columns[4])-1;
			if (columns[6] == "+")
			{
				gtfs.back().strand = 0;
			}
			else
			{
				gtfs.back().strand = 1;
			}
			vector<string> info;
			split(columns[8], ' ', info);
			for (int i = 0; i < info.size(); i += 2)
			{
				if (info[i] == "gene_id")
				{
					gtfs.back().gene_id = info[i + 1].substr(1, info[i + 1].size() - 3);
				}
				if (info[i] == "transcript_id")
				{
					gtfs.back().transcript_id = info[i + 1].substr(1, info[i + 1].size() - 3);
					//if(gtfs.back().transcript_id=="ENSDART00000143547")
					//{
					//	cout<<"\nline: "<<columns[8];
					//	getchar();
					//}
				}
				if (info[i] == "gene_biotype")
				{
					gtfs.back().gene_biotype = info[i + 1].substr(1, info[i + 1].size() - 3);
				}
				if (info[i] == "exon_number")
				{
					gtfs.back().exon_number = stoi(info[i + 1].substr(1, info[i + 1].size() - 3));
				}
				if(info[i] == "ccds_id")
				{
					gtfs.back().ccds_id = info[i + 1].substr(1, info[i + 1].size() - 3);
				}
			}
		}
	}
}

void read_gtf(vector<GTF> &gtfs, string filename, map<string,int> &chr_labels, bool includes_chr_prefix, int threads)
{
    // Read the whole file
    ifstream file(filename);
    stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    // Split the content by lines
    //string content = buffer.str();
    vector<string> lines;
    split(buffer.str(), '\n', lines);

    int num_lines = lines.size();

    // Process lines in parallel
    #pragma omp parallel num_threads(threads)
    {
        #pragma omp for schedule(dynamic) nowait
        for (int j = 0; j < num_lines; j++)
        {
            string &line = lines[j];
            if (line[0] != '#')
            {
                vector<string> columns;
                split(line, '\t', columns);
                if(!includes_chr_prefix)
                {
                    columns[0] = "chr"+columns[0];
                }
                if(!chr_labels.count(columns[0]))
                {
                    continue;
                }
                
                GTF newGtf;
                newGtf.chr = chr_labels[columns[0]];
                newGtf.contig = columns[0];
                newGtf.source = columns[1];
                newGtf.annotation_type = columns[2];
                newGtf.start = stoi(columns[3])-1;
                newGtf.end = stoi(columns[4])-1;
                newGtf.strand = (columns[6] == "+") ? 0 : 1;

                vector<string> info;
                split(columns[8], ' ', info);

                for (int i = 0; i < info.size(); i += 2)
                {
                    string key = info[i];
                    string value = info[i + 1].substr(1, info[i + 1].size() - 3);
                    
                    if (key == "gene_id") { newGtf.gene_id = value; }
                    if (key == "transcript_id") { newGtf.transcript_id = value; }
                    if (key == "gene_biotype") { newGtf.gene_biotype = value; }
                    if (key == "exon_number") { newGtf.exon_number = stoi(value); }
                    if (key == "ccds_id") { newGtf.ccds_id = value; }
                }

                #pragma omp critical
                {
                    gtfs.push_back(newGtf);
                }
            }
        }
    }
}

void construct_transcripts_original(vector<Transcript> &transcripts, const vector<GTF> &annotations)
{
	map<string, int> model_index;
	for (int i = 0; i < annotations.size(); i++)
	{
		if (annotations[i].annotation_type == "exon" && !model_index.count(annotations[i].transcript_id) && annotations[i].transcript_id!="")
		{
			model_index[annotations[i].transcript_id] = transcripts.size();
			transcripts.push_back(Transcript());
			transcripts.back().chr = annotations[i].chr;
			transcripts.back().strand = annotations[i].strand;
			transcripts.back().transcript_id = annotations[i].transcript_id;
		}
	}
	for (int i = 0; i < annotations.size(); i++)
	{
		if (annotations[i].annotation_type == "exon" && annotations[i].transcript_id!="")
		{
			int index = model_index.at(annotations[i].transcript_id);
			transcripts[index].exons.push_back(Exon());
			transcripts[index].exons.back().start = annotations[i].start;// -1;
			transcripts[index].exons.back().end = annotations[i].end;// -1;
			if (transcripts[index].exons.back().start < 0)
			{
				cout << "\ntranscript error: " << transcripts[index].exons.back().start << " " << i;
				getchar();
			}
		}
	}
}
void construct_transcripts(vector<Transcript> &transcripts, const vector<GTF> &annotations, int threads)
{
    map<string, int> model_index;
    // We need to synchronize threads when accessing shared resources
    omp_lock_t lock;
    omp_init_lock(&lock);

    #pragma omp parallel for num_threads(threads) shared(transcripts, model_index) schedule(dynamic)
    for (int i = 0; i < annotations.size(); ++i)
    {
        const GTF& annotation = annotations[i];
        if (annotation.annotation_type == "exon" && annotation.transcript_id != "")
        {
            omp_set_lock(&lock);
            // If the transcript_id does not exist in the map, create a new Transcript
            if (!model_index.count(annotation.transcript_id))
            {
                model_index[annotation.transcript_id] = transcripts.size();
                transcripts.push_back(Transcript());
                Transcript& transcript = transcripts.back();
                transcript.chr = annotation.chr;
                transcript.strand = annotation.strand;
                transcript.transcript_id = annotation.transcript_id;
            }

            // Now we know that the transcript_id exists in the map. Add a new Exon to it.
            int index = model_index.at(annotation.transcript_id);
            transcripts[index].exons.push_back(Exon());
            Exon& exon = transcripts[index].exons.back();
            exon.start = annotation.start;
            exon.end = annotation.end;
            omp_unset_lock(&lock);

        }
    }

    omp_destroy_lock(&lock);
}


void read_genome_speedy(vector<string> &genome, map<string, int> &chr_labels, string filename)
{
	
	ifstream file(filename);
	if(!file) {
		cerr << "Could not open " << filename << std::endl;
		exit(0);
	}

	string chr_str;
	string line;

	int chr_index = -1;
	while (getline(file, line))
	{
		if (line.substr(0,1) == ">")
		{

			vector<string> columns;
			split(line.substr(1), ' ', columns);
			chr_str = columns[0];

			chr_index++;
			chr_labels[chr_str] = chr_index;
			genome.push_back("");
		}
		else if (chr_index!=-1)
		{
			genome.back()+=line;
		}
	}
	
}



void read_genome_assign_reads(map<string, int> &chr_labels, string filename)
{
	ifstream file(filename);
	string chr_str;
	string line;
	int chr_index = -1;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{

			vector<string> columns;
			split(line.substr(1), ' ', columns);
			chr_str = columns[0];

			chr_index++;
			chr_labels[chr_str] = chr_index;
		}

	}
}

void get_transcript_seq_speedy(vector<Transcript> &transcripts, vector<string> &genome, int threads)
{
    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (int i = 0; i < transcripts.size(); i++)
    {
        Transcript& currTranscript = transcripts[i];  // Reference to current transcript

        for (int j = 0; j < currTranscript.exons.size(); j++)
        {
            Exon& currExon = currTranscript.exons[j];  // Reference to current exon

            currExon.seq = genome[currTranscript.chr].substr(currExon.start, currExon.end - currExon.start + 1);
            if (currTranscript.strand == 1) {
                //reverse(currExon.seq.begin(), currExon.seq.end());
				reverse_complement(currExon.seq);
            }
        }
    }
}


void get_orfs_genome(vector<CandidateORF> &orfs, vector<string> &genome, int threads, map<string,int> &chr_labels){
	
	//Make the reverse chr_labels
	map<int,string> rev_chr_labels;
	for(auto it=chr_labels.begin();it!=chr_labels.end();it++)
	{
		rev_chr_labels[it->second]=it->first;
	}
	

	//Holds both forward and reverse
	vector<string> double_genome;
	for(string sequence : genome){
		double_genome.emplace_back(sequence);
		reverse_complement(sequence);
		double_genome.emplace_back(sequence);
	}

	//Create the double genome, and use onwards
	
	
    threads = std::min(threads, static_cast<int>(double_genome.size()));
	cout <<"\nNum threads: " << to_string(omp_get_max_threads());

    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (int chr = 0; chr < double_genome.size(); chr++){
		map<int, int> prev_stops;

		
		//hi
        string& chromosome_seq = double_genome[chr];
		int seq_size = chromosome_seq.size();
		int strand = chr%2;
		
        for(int i=0; i<chromosome_seq.size()-2; i++){
            if(chromosome_seq.substr(i,3)=="ATG"){
                for(int j=i+3; j<chromosome_seq.size()-2; j+=3){
                    string codon = chromosome_seq.substr(j,3);
                    if(codon=="TAA" || codon=="TAG" || codon=="TGA"){
						
						//previous code had a length threshold of 7 MAKE THIS OPTIONAL
						if(j+2-i<7){
							break;
						}
						
                        CandidateORF orf = CandidateORF();
                        orf.chr = chr/2;

                        orf.orf_length = j+3-i;
						Exon fake_exon = Exon();
						
						if(strand == 0){
							fake_exon.start = i;
							fake_exon.end = j+2;
							orf.start_codon_pos = i;
							orf.stop_codon_pos = j+2;	
						}
						else if (strand == 1){
							//fake_exon.start = j+2;
							//fake_exon.end = i;
							//orf.start_codon_pos = j+2;
							//orf.stop_codon_pos = i;	
							
							fake_exon.start = seq_size - 3 - j;
							fake_exon.end = seq_size - 1 - i;
							orf.start_codon_pos = fake_exon.start;
							orf.stop_codon_pos = fake_exon.end;	
						}
						fake_exon.is_coding=true;
						fake_exon.seq = chromosome_seq.substr(i, j-i+3);
						orf.exons.emplace_back(fake_exon);
						orf.strand = strand;
						orf.chr_str = rev_chr_labels[chr/2];
						
						//Filter N sequences (MAKE THIS OPTIONAL)
						if(fake_exon.seq.find('N') != std::string::npos){
							break;
						}
						
						//Only longest sequences
						if(strand == 0 && prev_stops.find(orf.stop_codon_pos)!=prev_stops.end()){
							break;
						}
						else if(strand == 1 && prev_stops.find(orf.start_codon_pos)!=prev_stops.end()){
							break;
						}						
						
                        #pragma omp critical
                        {
                            orfs.emplace_back(orf);
							
							if(strand == 0){
								prev_stops[orf.stop_codon_pos]+=1;
							} else if(strand == 1){
								prev_stops[orf.start_codon_pos]+=1;
							}
                        }
                        break;
                    }
                }
            }
        }
    }
	
}


void get_orfs_speedy(vector<CandidateORF> &orfs, vector<Transcript> &transcripts, int threads)
{
	#pragma omp parallel for num_threads(threads) schedule(dynamic)
	for (int i = 0; i < transcripts.size(); i++)
	{
		Transcript& transcript = transcripts[i];
		sort(transcript.exons.begin(), transcript.exons.end());

		string transcript_seq = "";
		vector<int> exon_index;
		vector<int> genome_index;

		int total_exon_length = 0;
		for (Exon& exon : transcript.exons) {
			total_exon_length += exon.seq.length();
		}

		exon_index.reserve(total_exon_length);
		genome_index.reserve(total_exon_length);

		if(transcript.strand == 0){
			for(int j=0; j<transcript.exons.size(); j++){
				Exon& exon = transcript.exons[j];
				transcript_seq += exon.seq;
				for(int k=0; k<exon.seq.length(); k++){
					exon_index.emplace_back(j);
					genome_index.emplace_back(exon.start+k);
				}
			}
		}
		else if(transcript.strand == 1){
			for(int j=transcript.exons.size() - 1; j>=0; j--){
				Exon& exon = transcript.exons[j];
				transcript_seq += exon.seq;
				for(int k=exon.seq.length() - 1; k>=0; k--){
					exon_index.emplace_back(j);
					genome_index.emplace_back(exon.start+k);
				}
			}
		}



		for(int j=0; j<transcript_seq.size()-2; j++){
			if(transcript_seq.substr(j,3)=="ATG"){

				for(int k=j+3; k<transcript_seq.size()-2; k+=3){
										
					string codon = transcript_seq.substr(k,3);
					if(codon=="TAA" || codon=="TAG" || codon=="TGA"){
						
						
						CandidateORF orf = CandidateORF();
						
						orf.chr = transcript.chr;
						int start_exon = exon_index[j];
						int end_exon = exon_index[k+2];
						if (transcript.strand == 1)
						{
							start_exon= exon_index[k+2];
							end_exon = exon_index[j];
						}
						for (int w = start_exon; w <= end_exon; w++)
						{
							orf.exons.push_back(transcript.exons[w]);
							orf.exons.back().is_coding = true;
						}
						orf.exons.front().start = genome_index[j];// -1;
						orf.exons.back().end = genome_index[k + 2];// -1;
						orf.start_codon_pos = genome_index[j];// -1;
						orf.stop_codon_pos = genome_index[k + 2];// -1;
						if (transcript.strand == 1)
						{
							orf.exons.front().start = genome_index[k + 2]; 
							orf.exons.back().end = genome_index[j];
							orf.start_codon_pos = genome_index[k + 2];
							orf.stop_codon_pos = genome_index[j]; 

						}
						orf.strand = transcript.strand;
						orf.transcript_id = transcript.transcript_id;
						int orf_length = 0;
						for (int q = 0; q < orf.exons.size(); q++)
						{
							orf_length += 1+orf.exons[q].end - orf.exons[q].start;
						}
						orf.orf_length = orf_length;
						
						#pragma omp critical
						{
							orfs.emplace_back(orf);
						}
						break;
					}
				}
			}

			
		}

	}
	
}


void filter_duplicate_orfs(vector<CandidateORF> &orfs, int threads)
{
    set<string> orf_ids;
    vector<CandidateORF> filtered_orfs;
    map<int, map<int, int>> stop_locs;


    for (int i = 0; i < orfs.size(); i++)
    {
        CandidateORF &orf = orfs[i];
        if (orf.strand == 0)
        {
            if (!stop_locs.count(orf.chr) || !stop_locs.at(orf.chr).count(orf.stop_codon_pos) || orf.stop_codon_pos - orf.start_codon_pos > stop_locs.at(orf.chr).at(orf.stop_codon_pos))
            {

                stop_locs[orf.chr][orf.stop_codon_pos] = orf.stop_codon_pos - orf.start_codon_pos;
			   
            }
        }
        else
        {
            if (!stop_locs.count(orf.chr) || !stop_locs.at(orf.chr).count(orf.start_codon_pos) || orf.stop_codon_pos - orf.start_codon_pos > stop_locs.at(orf.chr).at(orf.start_codon_pos))
            {
				
					stop_locs[orf.chr][orf.start_codon_pos] = orf.stop_codon_pos - orf.start_codon_pos;
				
            }
        }
    }

    #pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (int i = 0; i < orfs.size(); i++)
    {
        CandidateORF &orf = orfs[i];
        string orf_id = to_string(orf.chr) + "_" + to_string(orf.start_codon_pos) + "_" + to_string(orf.stop_codon_pos);
        bool toAdd = false;
        if (orf.strand == 0 && stop_locs[orf.chr][orf.stop_codon_pos] == orf.stop_codon_pos - orf.start_codon_pos)
            toAdd = true;
        else if (orf.strand == 1 && stop_locs[orf.chr][orf.start_codon_pos] == orf.stop_codon_pos - orf.start_codon_pos)
            toAdd = true;

        if (toAdd)
        {
            #pragma omp critical
            {
                if (!orf_ids.count(orf_id))
                {
                    filtered_orfs.emplace_back(orf);
                    orf_ids.insert(orf_id);
                }
            }
        }
    }
    orfs = filtered_orfs;
}
void get_gene_seq_speedy(vector<CandidateORF> &genes, vector<string> &genome, int threads)
{
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < genes.size(); i++)
    {
        CandidateORF &gene = genes[i];
		string& genome_chr = genome[gene.chr];
        for (int j = 0; j < gene.exons.size(); j++)
        {
				
            Exon &exon = gene.exons[j];
			exon.seq = genome_chr.substr(exon.start, exon.end - exon.start + 1);

        }
    }
}





void expand_candidate_orfs(vector<CandidateORF> &gene_models, int threads)
{

	#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < gene_models.size(); i++)
	{
		CandidateORF &gene_model = gene_models[i];

		if ((gene_model.strand == 0 || gene_model.strand == 1) && gene_model.exons.size() > 0 && gene_model.stop_codon_pos != -1)
		{



			int cur_frame = 0;
			int start = (gene_model.strand == 0) ? 0 : gene_model.exons.size()-1;
			int end = (gene_model.strand == 0) ? gene_model.exons.size() : -1;
			int step = (gene_model.strand == 0) ? 1 : -1;

			for (int j = start; j != end; j += step)
			{
				Exon &exon = gene_model.exons[j];

				if (exon.end > gene_model.start_codon_pos && exon.start < gene_model.stop_codon_pos)
				{
					exon.is_coding = true;

					exon.start = max(exon.start, gene_model.start_codon_pos);
					exon.end = min(exon.end, gene_model.stop_codon_pos);
					
				}
			}
		}
	}
}

void expand_gene_models(vector<GeneModel> &gene_models, int threads)
{

	#pragma omp parallel for num_threads(threads)
	for (int i = 0; i < gene_models.size(); i++)
	{
		GeneModel &gene_model = gene_models[i];

		if ((gene_model.strand == 0 || gene_model.strand == 1) && gene_model.exons.size() > 0 && gene_model.stop_codon_pos != -1)
		{


			gene_model.frame = vector<int>(1 + gene_model.stop_codon_pos - gene_model.start_codon_pos, 4);

			int cur_frame = 0;
			int start = (gene_model.strand == 0) ? 0 : gene_model.exons.size()-1;
			int end = (gene_model.strand == 0) ? gene_model.exons.size() : -1;
			int step = (gene_model.strand == 0) ? 1 : -1;

			for (int j = start; j != end; j += step)
			{
				Exon &exon = gene_model.exons[j];

				if (exon.end > gene_model.start_codon_pos && exon.start < gene_model.stop_codon_pos)
				{
					exon.is_coding = true;

					exon.start = max(exon.start, gene_model.start_codon_pos);
					exon.end = min(exon.end, gene_model.stop_codon_pos);
					
					for (int k = 0; k <= exon.end - exon.start; k++)
					{
						int position = exon.start - gene_model.start_codon_pos + k;

						if (gene_model.strand == 1)
						{
							position = gene_model.stop_codon_pos - exon.end + k;
						}


						gene_model.frame[position] = cur_frame;

						cur_frame++;
						if (cur_frame == 3)
						{
							cur_frame = 0;
						}
					}
				}
			}
		}
	}
}


void print_genes(const vector<vector<CandidateORF>> &gene_models,string filename, int strand, int& orf_index)
{
	map<int, char> rev_nucmap = { { 0,'A' },{ 1,'C' },{ 2 ,'G'},{ 3,'T' },{ 4,'N' },{ 5,'-' },{6,'#'} };

    std::ofstream file(filename, std::ios::app); // Open the file in append mode
	for (int g = 0; g < gene_models.size(); g++)
	{
		for (int i = 0; i < gene_models[g].size(); i++)
		{
			if(gene_models[g][i].strand!=strand){
				continue;
			}
			file << "\n" << orf_index << " " << gene_models[g][i].species<<" "<< gene_models[g][i].transcript_id << " " << gene_models[g][i].gene_id << " " << gene_models[g][i].chr << " " << gene_models[g][i].strand << " " << gene_models[g][i].start_codon_pos << " " << gene_models[g][i].stop_codon_pos << " ";
			for (int j = 0; j < gene_models[g][i].exons.size(); j++)
			{
				file << gene_models[g][i].exons[j].start << ",";
			}
			file << " ";
			for (int j = 0; j < gene_models[g][i].exons.size(); j++)
			{
				file << gene_models[g][i].exons[j].end << ",";
			}
			file << " ";
			for (int j = 0; j < gene_models[g][i].exons.size(); j++)
			{
				file << gene_models[g][i].exons[j].is_coding << ",";
			}
			/*for (int j = 0; j < gene_models[g][i].read_count.size(); j++)
			{
				file << " " << gene_models[g][i].read_count[j];
			}*/
			
			file << " " << gene_models[g][i].num_agree << " " << gene_models[g][i].pos_agree<<" "<<gene_models[g][i].lncrna<<" "<<gene_models[g][i].processed_pseudogene<<" "<<gene_models[g][i].unprocessed_pseudogene<<" "<<gene_models[g][i].protein_coding<<" "<<gene_models[g][i].biotype_other<<" "<<gene_models[g][i].intersect_gene<<" "<<gene_models[g][i].orf_length<<" "<<gene_models[g][i].age_class<<" "<<gene_models[g][i].secondary_gene_id<<" "<<gene_models[g][i].exons.size()<<" "<<gene_models[g][i].intersect_five_utr<<" "<<gene_models[g][i].intersect_three_utr<<" "<<gene_models[g][i].uorf_id<<" "<<gene_models[g][i].dorf_id<<" "<<gene_models[g][i].antisense_gene<<" "<<gene_models[g][i].CDS_intersect<<" "<<gene_models[g][i].CDS_biotype<<" "<<gene_models[g][i].orf_annotation_type<<" "<<gene_models[g][i].orf_gene_biotype<<" "<<gene_models[g][i].chr_str;
			
			file << " ";
			for(Exon exon: gene_models[g][i].exons){
				file << exon.seq;
			}
			orf_index++;
		}
	}
}



void find_intersect_ann(vector<CandidateORF> &orfs, vector<GTF> &anns, map<string,CDS> &cds, int threads)
{
	int canonical_count =0;
	int ann_start_index = 0;
	int splice_count=0;
	for (int i = 0; i < orfs.size(); i++)
	{
		CandidateORF &my_orf = orfs[i];
		my_orf.gene_id="";
		while (anns[ann_start_index].chr < my_orf.chr || my_orf.start_codon_pos-anns[ann_start_index].end>1000000)
		{
			ann_start_index++;
		}
		for (int j = ann_start_index; j < anns.size(); j++)
		{
			GTF& my_ann = anns[j];
			if (my_ann.chr > my_orf.chr || my_ann.start > my_orf.stop_codon_pos)
			{
				break;
			}
			//intersect on opposite strand
			if (my_orf.chr == my_ann.chr && my_orf.strand != my_ann.strand && (
				(my_orf.start_codon_pos <= my_ann.start && my_orf.stop_codon_pos >= my_ann.start) ||
				(my_orf.start_codon_pos <= my_ann.end && my_orf.stop_codon_pos >= my_ann.end) ||
				(my_orf.start_codon_pos >= my_ann.start && my_orf.stop_codon_pos <= my_ann.end) ||
				(my_orf.start_codon_pos <= my_ann.start && my_orf.stop_codon_pos >= my_ann.end)))
			{
				if (my_ann.annotation_type == "exon")
				{
					my_orf.antisense_gene = my_ann.gene_id;
				}

			}
			//intersect on same strand
			if (my_orf.chr == my_ann.chr && my_orf.strand == my_ann.strand && (
				(my_orf.start_codon_pos <= my_ann.start && my_orf.stop_codon_pos >= my_ann.start) ||
				(my_orf.start_codon_pos <= my_ann.end && my_orf.stop_codon_pos >= my_ann.end) ||
				(my_orf.start_codon_pos >= my_ann.start && my_orf.stop_codon_pos <= my_ann.end) ||
				(my_orf.start_codon_pos <= my_ann.start && my_orf.stop_codon_pos >= my_ann.end)))
			{
				if (my_ann.annotation_type == "exon")
				{
					if (my_ann.gene_biotype == "lncRNA")
					{
						my_orf.lncrna = true;
					}
					else if (my_ann.gene_biotype == "transcribed_unprocessed_pseudogene"||my_ann.gene_biotype=="unprocessed_pseudogene")
					{
						my_orf.unprocessed_pseudogene = true;
					}
					else if (my_ann.gene_biotype == "protein_coding")
					{
						my_orf.protein_coding = true;
					}
					else if (my_ann.gene_biotype == "processed_pseudogene")
					{
						my_orf.processed_pseudogene = true;
					}
					else
					{
						my_orf.biotype_other = my_ann.gene_biotype;
					}
					my_orf.intersect_gene = my_ann.gene_id;
				}
				else if (my_ann.annotation_type == "five_prime_utr")
				{
					my_orf.intersect_five_utr = my_ann.gene_id;
				}
				else if (my_ann.annotation_type == "three_prime_utr")
				{
					my_orf.intersect_three_utr = my_ann.gene_id;
				}
				else if (my_ann.annotation_type == "CDS")
				{
					my_orf.CDS_intersect = my_ann.gene_id;
					if(my_ann.spliced_gene){
						splice_count++;
					}
					my_orf.CDS_biotype = my_ann.gene_biotype;
					string cds_id=my_ann.ccds_id;
					if(cds_id=="")
					{
						cds_id=my_ann.transcript_id;
					}
                    if(my_orf.strand==0)
                    {
                        if(abs(cds[cds_id].start-my_orf.start_codon_pos)<=3 && abs(cds[cds_id].end-my_orf.stop_codon_pos)<=3)
                        {   
							canonical_count++;
                            my_orf.gene_id = my_ann.gene_id;

							if(my_orf.gene_id==""){
								my_orf.gene_id="canonical";
							}
                            my_orf.orf_annotation_type = my_ann.annotation_type;
                            my_orf.orf_gene_biotype = my_ann.gene_biotype;
                        }
                    }
                    if(my_orf.strand==1)
                    {
                        if(abs(cds[cds_id].start-my_orf.start_codon_pos)<=3 && abs(cds[cds_id].end-my_orf.stop_codon_pos)<=3)
                        {   
							canonical_count++;

                            my_orf.gene_id = my_ann.gene_id;
							if(my_orf.gene_id==""){
								my_orf.gene_id="canonical";
							}
                            my_orf.orf_annotation_type = my_ann.annotation_type;
                            my_orf.orf_gene_biotype = my_ann.gene_biotype;
                        }
                    }
				}
			}
		}
	}
	cout << "\nCanonical genes: " << to_string(canonical_count);
	cout << "\nSpliced gene overlaps: " << to_string(splice_count);

}



void assemble_cds(map<string,CDS> &cds,const vector<GTF> &gtfs)
{
	map<string,map<int,int>> cds_exons;
	int n = gtfs.size();

	for(int i=0;i<n;i++)
	{
		const GTF& curr_gtf = gtfs[i];
		
		if(curr_gtf.annotation_type=="CDS")
		{
			string cds_id=curr_gtf.ccds_id;
			if(cds_id=="")
			{
				cds_id=curr_gtf.transcript_id;
			}
			int start = curr_gtf.start;
			int end = curr_gtf.end;
			int exon_number = curr_gtf.exon_number;

			CDS& curr_cds = cds[cds_id];
			curr_cds.exon_starts[exon_number]=start;
			curr_cds.exon_stops[exon_number]=end;
			if(curr_cds.start==-1 || start<curr_cds.start)
			{
				curr_cds.start=start;
			}
			if(curr_cds.end==-1 || end>curr_cds.end)
			{
				curr_cds.end=end;
			}
			
		}
	}
}




void read_sam_file(string filename, map<string, int> &chr_labels, int threads, vector<string> & sam_lines)
{
	string real_filename = filename;
	//Convert bam to sam, use sam

	std::string command_1 = "samtools view -@ " + std::to_string(threads) + " -h -o " + filename.substr(0, filename.length()-4) + ".sam " + filename;
	const char* command = command_1.c_str();
	if (filename.substr(filename.length()-4) == ".bam"){
			int a = system(command);
			filename = filename.substr(0, filename.length()-4) + ".sam";
	}
	

	ifstream file(filename);
	string line;
	int index =0;

	while (getline(file, line))
	{
		sam_lines.emplace_back(line);
	}

	//Convert bam to sam, use sam
	if (real_filename.substr(real_filename.length()-4) == ".bam"){
		system(("rm " + filename).c_str());
	}

}

void process_sam_lines(vector<Read> &reads, string filename, map<string, int> &chr_labels, vector<string> &sam_lines, int threads)
{
    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < sam_lines.size(); i++)
    {
        string line = sam_lines[i];
        if (line[0] != '@')
        {
            vector<string> cols;
            split(line, '\t', cols);
            string chr_str = cols[2];
            //if (chr_labels.count(chr_str))
            //{
                Read my_read = Read();
                my_read.chr = chr_labels[chr_str];
                my_read.start = stoi(cols[3]);

                int flag = stoi(cols[1]);
                my_read.flag = flag;
                    my_read.length = cols[9].size();

                if (flag & 16)
                {
                    my_read.strand = 1;
                    my_read.start += my_read.length-1;
                }
                else
                {
                    my_read.strand = 0;
                }

                #pragma omp critical
				{
					reads[i] = my_read;
				}
            //}
        }
        sam_lines[i] = "";
    }
}


/*
void read_genes_yeast(vector<GeneModel> & gene_models,string filename)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, ' ', column_data);
		GeneModel my_gene_model = GeneModel();
		//my_gene_model.species = column_data[1];
		//my_gene_model.transcript_id = column_data[2];
		if (column_data[14] == "Verified" || column_data[14] == "Uncharacterized" || column_data[14] == "transposable_element_gene"){
			my_gene_model.gene_id = column_data[14];
		} else{
				my_gene_model.gene_id = "";
		}
			my_gene_model.chr = stoi(column_data[2]);
			my_gene_model.strand = stoi(column_data[5]);
			my_gene_model.start_codon_pos = stoi(column_data[3]);
			my_gene_model.stop_codon_pos = stoi(column_data[4]);



			Exon my_exon = Exon();
			my_exon.start = my_gene_model.start_codon_pos;
			my_exon.end = my_gene_model.stop_codon_pos;
			//my_exon.is_coding = stoi(exon_coding[i]);
			my_gene_model.exons.emplace_back(my_exon);

			
			for (int i = 0; i < 3; i++)
			{//lncrna processed_pseudogene unprocessed_pseudogene protein_coding biotype_other intersect_gene_id

				my_gene_model.read_count[i] = stoi(column_data[11 + i]);
			}
			
			my_gene_model.num_agree = stoi(column_data[14]);
			my_gene_model.pos_agree = stoi(column_data[15]);
			my_gene_model.lncrna = stoi(column_data[16]);
			my_gene_model.processed_pseudogene = stoi(column_data[17]);
			my_gene_model.unprocessed_pseudogene = stoi(column_data[18]);
			my_gene_model.protein_coding = stoi(column_data[19]);
			my_gene_model.biotype_other = column_data[20];
			my_gene_model.intersect_gene = column_data[21];
			my_gene_model.orf_length = stoi(column_data[22]);
			
			gene_models.emplace_back(my_gene_model);

	}
}
*/
void read_genes(vector<GeneModel> & gene_models,string filename, bool annotated_only)
{
	ifstream file(filename);
	string line;
	getline(file, line);
	while (getline(file, line))
	{
		vector<string> column_data;
		split(line, ' ', column_data);
		
		//If the gene is noncanonical and we only want annotated, skip it
		if(annotated_only && column_data[3] == ""){
			continue;
		}
		GeneModel my_gene_model = GeneModel();
		//my_gene_model.species = column_data[1];
		//my_gene_model.transcript_id = column_data[2];
		my_gene_model.gene_id = column_data[3];
		my_gene_model.chr = stoi(column_data[4]);
		my_gene_model.strand = stoi(column_data[5]);
		my_gene_model.start_codon_pos = stoi(column_data[6]);
		my_gene_model.stop_codon_pos = stoi(column_data[7]);
		vector<string> exon_starts;
		split(column_data[8], ',', exon_starts);
		vector<string> exon_stops;
		split(column_data[9], ',', exon_stops);
		vector<string> exon_coding;
		split(column_data[10], ',', exon_coding);

		for (int i = 0; i < exon_starts.size(); i++)
		{
			Exon my_exon = Exon();
			my_exon.start = stoi(exon_starts[i]);
			my_exon.end = stoi(exon_stops[i]);
			my_exon.is_coding = stoi(exon_coding[i]);
			my_gene_model.exons.emplace_back(my_exon);

		}
		/*
		for (int i = 0; i < 3; i++)
		{//lncrna processed_pseudogene unprocessed_pseudogene protein_coding biotype_other intersect_gene_id

			my_gene_model.read_count[i] = stoi(column_data[11 + i]);
		}
		*/
		//my_gene_model.num_agree = stoi(column_data[14]);
		//my_gene_model.pos_agree = stoi(column_data[15]);
		//my_gene_model.lncrna = stoi(column_data[16]);
		//my_gene_model.processed_pseudogene = stoi(column_data[17]);
		//my_gene_model.unprocessed_pseudogene = stoi(column_data[18]);
		//my_gene_model.protein_coding = stoi(column_data[19]);
		//my_gene_model.biotype_other = column_data[20];
		//my_gene_model.intersect_gene = column_data[21];
		//my_gene_model.orf_length = stoi(column_data[22]);
		
		gene_models.emplace_back(my_gene_model);

	}
}
/*
void read_annotated_genes(vector<GeneModel>& genes, string filename)
{
	bool all_classes = true;
	const map<string, int> roman = { { "I",1 },{ "II",2 },{ "III",3 },{ "IV",4 },{ "V",5 },{ "VI",6 },{ "VII",7 },{ "VIII",8 },{ "IX",9 },{ "X",10 },{ "XI",11 },{ "XII",12 },{ "XIII",13 },{ "XIV",14 },{ "XV",15 },{ "XVI",16 },{ "Mito",17 } };
	ifstream file(filename);
	string line;
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',4 },{ 'n',4 } };
	bool read_pos = false;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{
			read_pos = false;
			string id = "";
			vector<string> column_data;
			split(line, ',', column_data);
			for (unsigned int i = 0; i < column_data.size(); i++)
			{
				vector<string> row_data;
				split(column_data[i], ' ', row_data);
				if (i == 0)
				{
					id = row_data[0].substr(1);
				}
				else if (row_data.size() > 1 && row_data[1] == "Chr")
				{
					string orf_class = "";
					if (line.find(", Verified ORF,") != string::npos)
					{
						orf_class = "verified";
					}
					else if (line.find(", Uncharacterized ORF,") != string::npos)
					{
						orf_class = "uncharacterized";
					}
					else if (line.find(", pseudogene,") != string::npos)
					{
						orf_class = "pseudogene";
					}
					else if (line.find(", telomere,") != string::npos)
					{
						orf_class = "telomere";
					}
					else if (line.find(", blocked_reading_frame,") != string::npos)
					{
						orf_class = "blocked";
					}
					else if ((line.find(", transposable_element_gene,") != string::npos))
					{
						orf_class = "te";
					}
					else if ((line.find(", LTR_retrotransposon,") != string::npos))
					{
						orf_class = "rte";
					}
					else if (line.find(", Dubious ORF,") != string::npos)// Dubious ORF
					{
						orf_class = "dubious";
					}

					if (all_classes || orf_class != "")
					{
						genes.push_back(GeneModel());
						genes.back().gene_id = orf_class;
						//genes.back().annotation = id;
						//genes.back().annotated = true;
						if (!roman.count(row_data[2]))
						{
							cout << "\nMISS: " << row_data[2];
							getchar();
						}
						genes.back().chr = roman.at(row_data[2]) - 1;
						int pos0 = stoi(row_data[4].substr(0, row_data[4].find('-'))) - 1;
						int pos1 = stoi(row_data[4].substr(row_data[4].find('-') + 1)) - 1;
						if (pos0 > pos1)
						{
							genes.back().strand = 1;
							genes.back().start_codon_pos = pos1;
							genes.back().stop_codon_pos = pos0;
						}
						else
						{
							genes.back().strand = 0;
							genes.back().start_codon_pos = pos0;
							genes.back().stop_codon_pos = pos1;
						}
						read_pos = true;
						Exon my_exon = Exon();
						my_exon.start = genes.back().start_codon_pos;
						my_exon.end = genes.back().stop_codon_pos;
						genes.back().exons.emplace_back(my_exon);
					}
				}
				else if (read_pos)
				{
					string pro_num = column_data[i].substr(0, column_data[i].find('-'));
					if (is_integer(pro_num))
					{
						int pos0 = stoi(pro_num) - 1;
						int pos1 = stoi(column_data[i].substr(column_data[i].find('-') + 1)) - 1;
						//genes.back().splice = 1;
						if (pos0 > pos1)
						{
							genes.back().stop_codon_pos = pos0;
							genes.back().exons[0].end = pos0;
						}
						else
						{
							genes.back().stop_codon_pos = pos1;
							genes.back().exons[0].end = pos1;

						}
					}
					else
					{
						break;
					}
				}
			}
		}
		
		else if (read_pos)
		{
			genes.back().seq_len += line.size();
			//cout << "\n " << to_string(genes.back().seq_len);
		}
	}
}
*/
void read_annotated_genes(vector<GeneModel>& genes, string filename)
{
	bool all_classes = true;
	const map<string, int> roman = { { "I",1 },{ "II",2 },{ "III",3 },{ "IV",4 },{ "V",5 },{ "VI",6 },{ "VII",7 },{ "VIII",8 },{ "IX",9 },{ "X",10 },{ "XI",11 },{ "XII",12 },{ "XIII",13 },{ "XIV",14 },{ "XV",15 },{ "XVI",16 },{ "Mito",17 } };
	ifstream file(filename);
	string line;
	map<char, int> nucmap = { { 'A',0 },{ 'a',0 },{ 'C',1 },{ 'c',1 },{ 'G',2 },{ 'g',2 },{ 'T',3 },{ 't',3 },{ 'N',4 },{ 'n',4 } };
	bool read_pos = false;
	while (getline(file, line))
	{
		if (line.substr(0, 1) == ">")
		{
			read_pos = false;
			string id = "";
			vector<string> column_data;
			split(line, ',', column_data);
			for (unsigned int i = 0; i < column_data.size(); i++)
			{
				vector<string> row_data;
				split(column_data[i], ' ', row_data);
				if (i == 0)
				{
					id = row_data[0].substr(1);
				}
				else if (row_data.size() > 1 && row_data[1] == "Chr")
				{
					string orf_class = "";
					if (line.find(", Verified ORF,") != string::npos)
					{
						orf_class = "verified";
					}
					else if (line.find(", Uncharacterized ORF,") != string::npos)
					{
						orf_class = "uncharacterized";
					}
					else if (line.find(", pseudogene,") != string::npos)
					{
						orf_class = "pseudogene";
					}
					else if (line.find(", telomere,") != string::npos)
					{
						orf_class = "telomere";
					}
					else if (line.find(", blocked_reading_frame,") != string::npos)
					{
						orf_class = "blocked";
					}
					else if ((line.find(", transposable_element_gene,") != string::npos))
					{
						orf_class = "te";
					}
					else if ((line.find(", LTR_retrotransposon,") != string::npos))
					{
						orf_class = "rte";
					}
					else if (line.find(", Dubious ORF,") != string::npos)// Dubious ORF
					{
						orf_class = "dubious";
					}

					if (all_classes || orf_class != "")
					{
						genes.push_back(GeneModel());
						genes.back().gene_id = orf_class;
						//genes.back().annotation = id;
						//genes.back().annotated = true;
						if (!roman.count(row_data[2]))
						{
							cout << "\nMISS: " << row_data[2];
							getchar();
						}
						genes.back().chr = roman.at(row_data[2]) - 1;
						int pos0 = stoi(row_data[4].substr(0, row_data[4].find('-'))) - 1;
						int pos1 = stoi(row_data[4].substr(row_data[4].find('-') + 1)) - 1;
						if (pos0 > pos1)
						{
							genes.back().strand = 1;
							genes.back().start_codon_pos = pos1;
							genes.back().stop_codon_pos = pos0;
						}
						else
						{
							genes.back().strand = 0;
							genes.back().start_codon_pos = pos0;
							genes.back().stop_codon_pos = pos1;
						}
						read_pos = true;
						Exon my_exon = Exon();
						my_exon.start = genes.back().start_codon_pos;
						my_exon.end = genes.back().stop_codon_pos;
						genes.back().exons.emplace_back(my_exon);
					}
				}
				else if (read_pos)
				{
					string pro_num = column_data[i].substr(0, column_data[i].find('-'));
					if (is_integer(pro_num))
					{
						int pos0 = stoi(pro_num) - 1;
						int pos1 = stoi(column_data[i].substr(column_data[i].find('-') + 1)) - 1;
						//genes.back().splice = 1;
						if (pos0 > pos1)
						{
							genes.back().stop_codon_pos = pos0;
							genes.back().exons[0].end = pos0;
						}
						else
						{
							genes.back().stop_codon_pos = pos1;
							genes.back().exons[0].end = pos1;

						}
					}
					else
					{
						break;
					}
				}
			}
		}
		
		else if (read_pos)
		{
			//genes.back().seq_len += line.size();
			//cout << "\n " << to_string(genes.back().seq_len);
		}
	}
}
void map_reads_to_genome_old(vector<vector<map<int, int>>> &reads_map_f, vector<vector<map<int, int>>> &reads_map_r, const vector<Read> &reads, int chromosome_count, int threads)
{
    int count = 0;
    int rm_size = reads_map_f.size();
    reads_map_f = vector<vector<map<int, int>>>(40, vector<map<int, int>>(chromosome_count));
    reads_map_r = vector<vector<map<int, int>>>(40, vector<map<int, int>>(chromosome_count));

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < reads.size(); i++)
    {
        Read my_read = reads[i];
        if (my_read.length >= 0 && my_read.length < rm_size)
        {
            if (my_read.strand == 0)
            {
                #pragma omp critical
				{
					reads_map_f[my_read.length][my_read.chr][my_read.start] += 1;
				}
            }
            else if (my_read.strand == 1)
            {
                #pragma omp critical
				{				
					reads_map_r[my_read.length][my_read.chr][my_read.start] += 1;
				}
            }
        }
    }
}
void map_reads_to_genome(vector<vector<map<int, int>>> &reads_map_f, vector<vector<map<int, int>>> &reads_map_r, const vector<Read> &reads, int chromosome_count, int threads)
{
    int count = 0;
    int rm_size = reads_map_f.size();
    reads_map_f = vector<vector<map<int, int>>>(40, vector<map<int, int>>(chromosome_count));
    reads_map_r = vector<vector<map<int, int>>>(40, vector<map<int, int>>(chromosome_count));

    vector<vector<map<int, int>>> local_reads_map_f[threads];
    vector<vector<map<int, int>>> local_reads_map_r[threads];
    
    for (int i = 0; i < threads; ++i)
    {
        local_reads_map_f[i] = vector<vector<map<int, int>>>(40, vector<map<int, int>>(chromosome_count));
        local_reads_map_r[i] = vector<vector<map<int, int>>>(40, vector<map<int, int>>(chromosome_count));
    }

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < reads.size(); i++)
    {
        int tid = omp_get_thread_num();
        Read my_read = reads[i];
        if (my_read.length >= 0 && my_read.length < rm_size)
        {
            if (my_read.strand == 0)
            {
                local_reads_map_f[tid][my_read.length][my_read.chr][my_read.start] += 1;
            }
            else if (my_read.strand == 1)
            {
                local_reads_map_r[tid][my_read.length][my_read.chr][my_read.start] += 1;
            }
        }
    }

	// Combine results from all threads
	for (int t = 0; t < threads; ++t)
	{
		for (int len = 0; len < 40; ++len)
		{
			auto& len_map_f = reads_map_f[len];
			auto& len_map_r = reads_map_r[len];

			for (int chr = 0; chr < chromosome_count; ++chr)
			{
				auto& chr_map_f = len_map_f[chr];
				auto& chr_map_r = len_map_r[chr];

				#pragma omp parallel sections num_threads(2)
				{
					#pragma omp section
					{
						for (auto& item : local_reads_map_f[t][len][chr])
						{
							chr_map_f[item.first] += item.second;
						}
					}
					#pragma omp section
					{
						for (auto& item : local_reads_map_r[t][len][chr])
						{
							chr_map_r[item.first] += item.second;
						}
					}
				}
			}
		}
	}

}

int find_p_site(vector<GeneModel> &orfs, vector<map<int, int>> &reads_map)
{
	int nearby_len = 40;
	int displacement = 20;
	int total_reads = 0;

	vector<int> nearby(nearby_len);
	vector<int> frame(3);
	int best_frame = -1;
	int max_frame = 0;
	int max_pos_reads = 0;
	int predicted_p_site = -1000;


	for (int i = 0; i < orfs.size(); i++)
	{
		GeneModel &my_orf = orfs[i];
		if (my_orf.strand==0)
		{
			map<int,int> &my_map = reads_map[my_orf.chr];
			vector<int> nearby_frames(3);

			for (int j = 0; j < nearby_len; j++)
			{
				int pos = my_orf.start_codon_pos - displacement + j;
				auto it = my_map.find(pos);
				if (it != my_map.end())
				{
					// If the element is found
					int value = it->second;

						nearby[j] += value;
						nearby_frames[j%3] += value;
						total_reads += value;
					
				}
			}
		}
	}


	for (int i = 0; i < nearby.size(); i++)
	{
		frame[i % 3] += nearby[i];
		if (nearby[i] > max_pos_reads)
		{
			max_pos_reads = nearby[i];
		}
	}


	for (int i = 0; i < frame.size(); i++)
	{
		if (frame[i] > max_frame)
		{
			best_frame = i;
			max_frame = frame[i];
		}
	}

	for (int i = 0; i < nearby.size(); i++)
	{
		//if (i % 3 == best_frame && nearby[i] > total_reads*.05)
		if (i%3 == best_frame && nearby[i] > max_pos_reads * .35)		
		{
			predicted_p_site = i - displacement;
			break;
		}
	}

	return predicted_p_site;
}

int find_p_site_rc(vector<GeneModel> &orfs, vector<map<int, int>> &reads_map)
{
	int nearby_len = 40;
	int displacement = 20;
	int total_reads = 0;

	vector<int> nearby(nearby_len);
	vector<int> frame(3);
	int best_frame = -1;
	int max_frame = 0;
	int max_pos_reads = 0;
	int predicted_p_site = -1000;


	for (int i = 0; i < orfs.size(); i++)
	{
		GeneModel &my_orf = orfs[i];
		if (my_orf.strand==1)
		{
			map<int,int> &my_map = reads_map[my_orf.chr];
			vector<int> nearby_frames(3);

			for (int j = 0; j < nearby_len; j++)
			{
				int pos = my_orf.stop_codon_pos - displacement + j;
				auto it = my_map.find(pos);
				if (it != my_map.end())
				{
					int value = it->second;

					nearby[j] += value;
					nearby_frames[j%3] += value;
					total_reads += value;
				}
			}
		}
	}


	for (int i = nearby.size()-1; i >= 0; i--)
	{
		frame[i % 3] += nearby[i];
		if (nearby[i] > max_pos_reads)
		{
			max_pos_reads = nearby[i];
		}
	}


	for (int i = 0; i < frame.size(); i++)
	{
		if (frame[i] > max_frame)
		{
			best_frame = i;
			max_frame = frame[i];
		}
	}

	for (int i = nearby.size()-1; i >= 0; i--)
	{
		//if (i % 3 == best_frame && nearby[i] > total_reads*.05)
		if (i%3 == best_frame && nearby[i] > max_pos_reads * .35)
		{
			predicted_p_site = i - displacement;
			break;
		}
	}

	return predicted_p_site;
}




void print_gene_reads_new(vector<GeneModel> &genes,string filename)
{
	ofstream file(filename);
	file << "gene_id reads";
	string temp = "";
	int lines = 0;
	cout << "\nGENES SIZE: " << to_string(genes.size());
	
	for (int i = 0; i < genes.size(); i++)
	{
		int printed_line = 0;
		for (int j = 0; j < genes[i].pos_read_count.size(); j++)
		{
			int num = genes[i].pos_read_count[j];
			if (num !=0){
				if (printed_line==0){
					temp+= "\n" + to_string(i) + " ";
					printed_line = 1;

				}
				//file << j << "-" << num << ",";
				temp += to_string(j) + "-" + to_string(num) + "-" + to_string(genes[i].frame[j]) + ",";
				lines++;
				if (lines == 100){
					file << temp;
					lines = 0;
					temp = "";
				}
			}
		}
	}
	file << temp;
}
void filter_orfs(vector<CandidateORF> &gene_models, int threads)
{
	vector<CandidateORF> filtered_gene_models;
	filtered_gene_models.reserve(gene_models.size());


	#pragma omp parallel for num_threads(threads)
	for(int i = 0; i < gene_models.size(); i++)
	{
		CandidateORF &my_gene_model = gene_models[i];
		if (1 + my_gene_model.stop_codon_pos - my_gene_model.start_codon_pos >= 0 && my_gene_model.stop_codon_pos - my_gene_model.start_codon_pos <= 10000000)
		{
			#pragma omp critical
			{
				filtered_gene_models.emplace_back(my_gene_model);
			}
		}
	}

	gene_models = move(filtered_gene_models);
}

void read_riboseq_studies(vector<Study> &studies, string filename)
{
	ifstream file(filename);
	string line;
	while (getline(file, line))
	{
		if(line.size()>0)
		{
			vector<string> columns;
			split(line, ' ', columns);
			studies.push_back(Study());
			studies.back().path = columns[0];
		}        
		//studies.back().srr = columns[1];
		//studies.back().preprocessed = stoi(columns[3]);
		//if (!studies.back().preprocessed)
		//{
		//	studies.back().linker = columns[4];
		//}
	}
}




void get_orf_chr_str(vector<CandidateORF> &orfs, map<string,int> &chr_labels, int threads)
{
	map<int,string> rev_chr_labels;
	for(auto it=chr_labels.begin();it!=chr_labels.end();it++)
	{
		rev_chr_labels[it->second]=it->first;
	}

    #pragma omp parallel for num_threads(threads) schedule(dynamic)
	for(int i=0;i<orfs.size();i++)
	{
		orfs[i].chr_str=rev_chr_labels.at(orfs[i].chr);
	}
}




bool quality_control(vector<GeneModel> &orfs, vector<map<int, int>> &reads_map_f, vector<map<int, int>> &reads_map_r, int p_site_f, int p_site_r, int threads, bool qc_positions) {
    vector<int> ann_frames(3, 0);

    #pragma omp parallel for num_threads(threads)
    for (int i = 0; i < orfs.size(); i++) {
        GeneModel &my_orf = orfs[i];
        vector<int> local_ann_frames(3, 0);

        if (my_orf.strand == 0) {
			map<int,int> &my_map = reads_map_f[my_orf.chr];
			auto my_map_end = my_map.end();
			
			for (int j = my_orf.start_codon_pos; j <= my_orf.stop_codon_pos; j++) {
				int frame_val = my_orf.frame[j - my_orf.start_codon_pos];
				if(frame_val > 2){
					continue;
				}

				int j_plus_p = j+p_site_f;

				if (my_map.find(j_plus_p)!= my_map_end) {
					if(!qc_positions){
						local_ann_frames[frame_val]+= my_map[j_plus_p];
					} else{
						local_ann_frames[frame_val]+= 1;
					}
				}
			}
        } else if (my_orf.strand == 1){
            map<int,int> &my_map = reads_map_r[my_orf.chr];
            auto my_map_end = my_map.end();

			for (int j = my_orf.stop_codon_pos; j >= my_orf.start_codon_pos; j--) {
				int frame_val = my_orf.frame[my_orf.stop_codon_pos - j];
				if(frame_val > 2){
					continue;
				}

				int j_plus_p = j+p_site_r;

				if (my_map.find(j_plus_p)!= my_map_end) {
					if(!qc_positions){
						local_ann_frames[frame_val]+= my_map[j_plus_p];
					} else{
						local_ann_frames[frame_val]+= 1;
					}
				}
			}
		}

        #pragma omp critical
        {
            for (int j=0; j<3; j++) {
                ann_frames[j] += local_ann_frames[j];
            }
        }
    }

    cout << "\nAnn frames: " << to_string(ann_frames[0]) << " " << to_string(ann_frames[1]) << " " << to_string(ann_frames[2]);
    
    //std::sort(ann_frames.begin(), ann_frames.end(), std::greater<int>());
    if(ann_frames[0]< 10000){
        return false;
    }

    //Check if frame0 is big enough
    if (ann_frames[0] < (ann_frames[1]) * 2.0 || ann_frames[0] < (ann_frames[2]) * 2.0){
        return false;
    }
    
    return true;
}

int findMaxIndex(int num1, int num2, int num3) {
    if (num1 > num2 && num1 > num3) {
        return 0;
    } else if (num2 > num1 && num2 > num3) {
        return 1;
    } else if (num3 > num1 && num3> num1){
        return 2;
    }
	
	return 3;
}
double log_nChoosek(unsigned n, unsigned k)
{
    if (k == 0) return 0;

    double result = log(n);
    for (int i = 2; i <= k; ++i) {
        result += log(n - i + 1);
        result -= log(i);
    }
    return result;
}

double binom_test(int k, int n, double p)
{
    double pval = 0;
    for (int i = k; i <= n; i++) // use 'i' instead of 'k' for the loop variable
    {
        double log_single_prob = log_nChoosek(n, i) + i * log(p) + (n - i) * log(1.0 - p);
        pval += exp(log_single_prob);
    }
    return pval;
}

void assign_reads_to_orfs(vector<GeneModel> &all_orfs, vector<map<int, int>> &reads_map_f, vector<map<int, int>> &reads_map_r, int threads) {
   
    // Seed with a fixed number
    default_random_engine engine(42);
	
	
	#pragma omp parallel for num_threads(threads)
    for (int i = 0; i < all_orfs.size(); i++) {
        GeneModel &my_orf = all_orfs[i];


		vector<int> exon_pos_reads (my_orf.stop_codon_pos - my_orf.start_codon_pos + 1);
		my_orf.pos_read_count = vector<int>(my_orf.stop_codon_pos - my_orf.start_codon_pos + 1);



		if (my_orf.strand == 0) {
			
			map<int,int> &my_map = reads_map_f[my_orf.chr];
			auto my_map_end = my_map.end();
			
			for (int j = my_orf.start_codon_pos; j <= my_orf.stop_codon_pos; j++) {
				if(my_orf.frame[j-my_orf.start_codon_pos] > 2){
					continue;
				}


				if (my_map.find(j)!= my_map_end) {
					exon_pos_reads[j-my_orf.start_codon_pos] = (my_map[j]);
					my_orf.pos_read_count[j-my_orf.start_codon_pos]=my_map[j];
				}

			}
		} else if(my_orf.strand==1) {
			map<int,int> &my_map = reads_map_r[my_orf.chr];
			auto my_map_end = my_map.end();
			
			for (int j = my_orf.stop_codon_pos; j >= my_orf.start_codon_pos; j--) {
				if(my_orf.frame[my_orf.stop_codon_pos-j]>2){
					continue;
				}


				if (my_map.find(j)!= my_map_end) {
					exon_pos_reads[my_orf.stop_codon_pos-j] = my_map[j];
					my_orf.pos_read_count[my_orf.stop_codon_pos-j]=my_map[j];

				}

			}
		}
        

		//Count triplet pattern
		for(int j=0; j<exon_pos_reads.size()-2; j+=3){
			//int index = findMaxIndex(exon_pos_reads[j], exon_pos_reads[j+1], exon_pos_reads[j+2]);
			int first = exon_pos_reads[j];
			int second = exon_pos_reads[j+1];
			int third = exon_pos_reads[j+2];

			//if(index!=3){ //indicates there was a max index
				//my_orf.max_codon[findMaxIndex(exon_pos_reads[j], exon_pos_reads[j+1], exon_pos_reads[j+2])]+=1; //Find max codon position

			if(first>second && first>third){
				my_orf.first_max++;
			}
			if(first+second+third!=0){
				my_orf.total_frames++;
			}
			//}
			my_orf.frame_reads[0]+=first;
			my_orf.frame_reads[1]+=second;
			my_orf.frame_reads[2]+=third;
		}
		
		//Find scrambled data
		shuffle(exon_pos_reads.begin(), exon_pos_reads.end(), engine);
		//Count triplet pattern
		for(int j=0; j<exon_pos_reads.size()-2; j+=3){
			//int index = findMaxIndex(exon_pos_reads[j], exon_pos_reads[j+1], exon_pos_reads[j+2]);
			int first = exon_pos_reads[j];
			int second = exon_pos_reads[j+1];
			int third = exon_pos_reads[j+2];

			//if(index!=3){ //indicates there was a max index
				//my_orf.max_codon[findMaxIndex(exon_pos_reads[j], exon_pos_reads[j+1], exon_pos_reads[j+2])]+=1; //Find max codon position
			if(first>second && first>third){
				my_orf.first_max_scrambled++;
			}
			if(first+second+third!=0){
				my_orf.total_frames_scrambled++;
			}
			//}
			my_orf.frame_reads_scrambled[0]+=first;
			my_orf.frame_reads_scrambled[1]+=second;
			my_orf.frame_reads_scrambled[2]+=third;
		}
		/*
		for(int j=0; j<exon_pos_reads.size()-2; j+=3){
			int index = findMaxIndex(exon_pos_reads[j], exon_pos_reads[j+1], exon_pos_reads[j+2]);
			if(index!=3){ //indicates there was a max index
				my_orf.max_codon_scrambled[findMaxIndex(exon_pos_reads[j], exon_pos_reads[j+1], exon_pos_reads[j+2])]+=1; //Find max codon position
			}		
		}*/

		/*
		int frame_sum = my_orf.max_codon[0] + my_orf.max_codon[1] + my_orf.max_codon[2];
		int scrambled_sum = my_orf.max_codon_scrambled[0] + my_orf.max_codon_scrambled[1] + my_orf.max_codon_scrambled[2];

		
		if (my_orf.max_codon[0]>0){
			my_orf.p_value = binom_test(my_orf.max_codon[0], frame_sum, 1.0/3.0);
		} else{
			my_orf.p_value = 1;
		}
		
		if(my_orf.max_codon_scrambled[0]>0){
			my_orf.p_value_scrambled = binom_test(my_orf.max_codon_scrambled[0], scrambled_sum, 1.0/3.0);
		} else{
			my_orf.p_value_scrambled = 1;
		}
		*/
		//Clear orf frames for memory efficiency
		//my_orf.frame = vector<int>();
    }
	
	
}




void print_all_passed_reads(const vector<map<int,int>>& passed_reads, string filename, int strand)
{
    ofstream file(filename);
    file << "chr strand pos count";

    for(size_t i = 0; i < passed_reads.size(); ++i)
    {
        const auto& current_map = passed_reads[i];
        for(const auto& pair : current_map)
        {
            file << "\n" << i << " " << strand << " " << pair.first << " " << pair.second;
        }
    }
}



void accumulate_reads(int threads, 
                      map<string, int>& chr_labels, 
                      vector<map<int, int>>& reads_map_f, 
                      vector<map<int, int>>& reads_map_r, 
                      vector<map<int, int>>& all_passed_reads_f, 
                      vector<map<int, int>>& all_passed_reads_r,
					  int p_site_f,
					  int p_site_r)
{
    // Set the number of threads
    omp_set_num_threads(threads);

    for (int chr = 0; chr < chr_labels.size(); chr++)
    {
        auto& all_passed_reads_f_chr = all_passed_reads_f[chr];
        auto& all_passed_reads_r_chr = all_passed_reads_r[chr];

		auto& reads_map_f_chr = reads_map_f[chr];
		auto& reads_map_r_chr = reads_map_r[chr];

		#pragma omp parallel sections num_threads(2)
		{
			#pragma omp section
			{

				for (auto it = reads_map_f_chr.begin(); it != reads_map_f_chr.end(); ++it)
				{
					all_passed_reads_f_chr[it->first-p_site_f] += it->second;
				}
				
			}
			#pragma omp section
			{

				for (auto it = reads_map_r_chr.begin(); it != reads_map_r_chr.end(); ++it)
				{
					all_passed_reads_r_chr[it->first-p_site_r] += it->second;
				}
				
			}
		}
		
	}
}

void outputTracks(vector<map<int,int>>& passed_reads_f, vector<map<int,int>>& passed_reads_r, vector<GeneModel>& orfs, map<string, int>& chr_labels, string& output_dir){

    // Construct rev_chr_labels once
    map<int,string> rev_chr_labels;
    for(auto& kv : chr_labels)
    {
        rev_chr_labels[kv.second] = kv.first;
    }
    
    ofstream file_plus(output_dir + "riboseq_reads_plus.wig");
    ofstream file_minus(output_dir + "riboseq_reads_minus.wig");
    ofstream file(output_dir + "candidate_orfs.gff3");
    
    string buffer_plus, buffer_minus, buffer_file;
        
    for(auto& kv : chr_labels)
    {
        buffer_plus += "variableStep chrom=" + kv.first + "\n";
        for(auto& kv2 : passed_reads_f[kv.second])
        {
            buffer_plus += std::to_string(kv2.first) + " " + std::to_string(kv2.second) + "\n";
        }
        file_plus << buffer_plus;
        buffer_plus = ""; //clear the buffer

        buffer_minus += "variableStep chrom=" + kv.first + "\n";
        for(auto& kv2 : passed_reads_r[kv.second])
        {
            buffer_minus += std::to_string(kv2.first) + " " + std::to_string(kv2.second) + "\n";
        }
        file_minus << buffer_minus;
        buffer_minus = ""; //clear the buffer
    }
    
    for(int i=0;i<orfs.size();i++)
    {
        GeneModel& my_orf = orfs[i];
        for(int j=0;j<my_orf.exons.size();j++)
        {
            buffer_file += rev_chr_labels[my_orf.chr] + "\tiRibo\tCDS\t" + std::to_string(my_orf.exons[j].start) + "\t" + std::to_string(my_orf.exons[j].end) + "\t.\t";
            buffer_file += (my_orf.strand==0 ? "+" : "-");
            buffer_file += "\t0\tID=candidate_orf" + std::to_string(i) + "\n";
        }
        file << buffer_file;
        buffer_file = ""; //clear the buffer
    }   
}
void print_translation_calls(vector<GeneModel>& orfs, string filename){
	ofstream file(filename);
	file << "index frame0 frame_sum scrambled0 scrambled_sum reads0 reads1 reads2 scrambled_reads0 scrambled_reads1 scrambled_reads2";
	for(int i=0; i<orfs.size(); i++){
		GeneModel& my_orf = orfs[i];

		vector<int>& reads = my_orf.frame_reads;
		vector<int>& reads_scrambled = my_orf.frame_reads_scrambled;


		//Only print orfs with at least one read
		//if(reads[0]>0 || reads[1]>0 || reads[2]>0){
			file << "\n" << i << " " << my_orf.first_max << " " << my_orf.total_frames  << " " << my_orf.first_max_scrambled << " " << my_orf.total_frames_scrambled << " " << reads[0] << " " << reads[1] << " " << reads[2] << " " << reads_scrambled[0] << " " << reads_scrambled[1] << " " << reads_scrambled[2];// << " " << my_orf.p_value << " " << my_orf.p_value_scrambled;
		//}
	}
	
}
void create_dir(string output_dir){
	if (fs::exists(output_dir)) {
		std::cout << "Directory already exists." << std::endl;
	} else {
		// Attempt to create the directory
		if (fs::create_directory(output_dir)) {
			std::cout << "Directory created successfully." << std::endl;
		} else {
			std::cout << "Failed to create directory." << std::endl;
		}
	}	
}
bool stringToBool(const std::string& str) {
    // Convert the string to lowercase for case-insensitive comparison
    std::string lowercaseStr;
    lowercaseStr.reserve(str.size());
    for (char c : str) {
        lowercaseStr.push_back(tolower(c));
    }

    // Check for common representations of boolean values
    if (lowercaseStr == "true" || lowercaseStr == "1" || lowercaseStr == "yes" || lowercaseStr == "on") {
        return true;
    } else if (lowercaseStr == "false" || lowercaseStr == "0" || lowercaseStr == "no" || lowercaseStr == "off") {
        return false;
    }

    // If the string doesn't match any known representation, you can decide the default value or handle the error as needed.
    // In this example, the default value is set to false.
    return false;
}
int main(int argc, char *argv[])
{
	unordered_map<string, string> parseArgs = parseArguments(argc, argv);
	vector<string> args(argv + 1, argv + argc);

	if(args.size() == 0){
		cerr << "No arguments provided";
		exit(0);
	}
	
	string runMode = getArg(parseArgs, "RunMode", "", true);


if(runMode=="GetCandidateORFs")
{
    auto start = std::chrono::high_resolution_clock::now();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed;
    
	int threads = stoi(getArg(parseArgs, "Threads", "1", false));
	cout <<"\nNum threads: " << to_string(threads);
	string output_dir = getArg(parseArgs, "Output", "", false);
	if(output_dir.size()!=0){
		create_dir(output_dir);
		output_dir+="/";
	}
	
    cout << "RUNNING GCO\n";;
    string genome_path = getArg(parseArgs, "Genome", "", true);
    string genome_annotation_path = getArg(parseArgs, "Transcriptome", "", false);
    string canonical_gene_annotation_path = getArg(parseArgs, "Annotations", "", true);
    vector<string> genome;
    map<string, int> chr_labels;
    
	bool genomeOnly = (genome_annotation_path.size()==0);
	bool gff = (canonical_gene_annotation_path.size() >= 3 &&
            (canonical_gene_annotation_path.substr(canonical_gene_annotation_path.size() - 3) == "gff" ||
             (canonical_gene_annotation_path.size() >= 4 &&
              canonical_gene_annotation_path.substr(canonical_gene_annotation_path.size() - 4) == "gff3")));
			  
	cout << "\nreading genome at: " <<genome_path;
    start = std::chrono::high_resolution_clock::now();
    read_genome_speedy(genome, chr_labels, genome_path);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nread_genome took: " << elapsed.count() << " s";
    cout << "\ncontigs read: "<<genome.size();
    
	vector<CandidateORF> orfs;

    vector<GTF> annotations2; 
    start = std::chrono::high_resolution_clock::now();
	cout << "\n" << canonical_gene_annotation_path;
	
    //read_gtf_original(annotations2, args[3] , chr_labels, true);
	read_gff3(annotations2, canonical_gene_annotation_path, chr_labels);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nread_gtf for annotations2 took: " << elapsed.count() << " s";
    cout << "\nannotations read from file: " << canonical_gene_annotation_path[3] <<" "<< annotations2.size();
    start = std::chrono::high_resolution_clock::now();
    sort(annotations2.begin(),annotations2.end());
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nsort for annotations2 took: " << elapsed.count() << " s";
    
    map<string,CDS> cds;
    start = std::chrono::high_resolution_clock::now();
    assemble_cds(cds,annotations2);//, threads);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nassemble_cds took: " << elapsed.count() << " s";
    
	
	if(!genomeOnly){
		vector<GTF> annotations;
		cout << "\nreading annotations at: "<<genome_annotation_path;
		start = std::chrono::high_resolution_clock::now();
		//read_gtf(annotations, genome_annotation_path , chr_labels, true, threads);
		read_gtf_original(annotations, genome_annotation_path , chr_labels, true);

		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		cout << "\nread_gtf took: " << elapsed.count() << " s";
		cout << "\nannotations read: " << annotations.size();
		
		
		vector<Transcript> transcripts;
		start = std::chrono::high_resolution_clock::now();
		construct_transcripts_original(transcripts, annotations);
		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		cout << "\nconstruct_transcripts took: " << elapsed.count() << " s";
		cout << "\nconstruct transcripts: " << transcripts.size();
		
		annotations.clear();

		cout << "\nget transcript seq";
		start = std::chrono::high_resolution_clock::now();
		get_transcript_seq_speedy(transcripts, genome, threads);
		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		cout << "\nget_transcript_seq took: " << elapsed.count() << " s";


		start = std::chrono::high_resolution_clock::now();
		get_orfs_speedy(orfs, transcripts, threads);
		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		cout << "\nget_orfs took: " << elapsed.count() << " s";

		transcripts.clear();
	} else{
		start = std::chrono::high_resolution_clock::now();
		get_orfs_genome(orfs, genome, threads, chr_labels);
		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		cout << "\nget_orfs took: " << elapsed.count() << " s";
	}

    start = std::chrono::high_resolution_clock::now();
    //get_orf_chr_str(orfs,chr_labels, threads);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nget_orf_chr_str took: " << elapsed.count() << " s";
    cout << "\nget orfs: " << orfs.size();

    start = std::chrono::high_resolution_clock::now();
    //filter_duplicate_orfs(orfs, threads);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nfilter_duplicate_orfs took: " << elapsed.count() << " s";
    cout << "\norfs remaining after same-stop filtering: " << orfs.size();

    start = std::chrono::high_resolution_clock::now();
    //get_gene_seq_speedy(orfs, genome, threads);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nget_gene_seq took: " << elapsed.count() << " s";
    cout << "\nget gene seqs";
	
	genome.clear();
    

    start = std::chrono::high_resolution_clock::now();
    sort(orfs.begin(), orfs.end());

    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nsort took: " << elapsed.count() << " s";
    cout << "\nsort ORFs";

    start = std::chrono::high_resolution_clock::now();
    //filter_orfs(orfs, threads); DO THIS IN GET ORFS
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nfilter_orfs took: " << elapsed.count() << " s";
    
    start = std::chrono::high_resolution_clock::now();
	//expand_candidate_orfs(orfs, threads); DO THIS IN GET_ORFS
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nexpand_gene_models took: " << elapsed.count() << " s";
    cout << "\nexpand gene models: "<<orfs.size();;

    start = std::chrono::high_resolution_clock::now();
    find_intersect_ann(orfs, annotations2, cds, threads);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nfind_intersect_ann took: " << elapsed.count() << " s";

	/*for(GTF& ann :annotations2){
		if(ann.gene_id!=""){
			cout << "\n" << ann.gene_id;
		}
	}*/
	annotations2.clear();
	cds.clear();
	
    vector<vector<CandidateORF>> gene_models_shell;
    gene_models_shell.emplace_back(orfs);
    start = std::chrono::high_resolution_clock::now();
    ofstream file(output_dir + "candidate_orfs"); // Open the file in append mode
	file << "id species transcript_id gene_id chr strand start stop exon_starts exon_stops exon_coding num_agree pos_agree lncrna processed_pseudogene unprocessed_pseudogene protein_coding biotype_other intersect_gene_id orf_length age_class secondary_gene_id exon_count intersect_five_utr intersect_three_utr uorf_id dorf_id antisense_gene CDS_intersect CDS_biotype orf_annotation_type orf_gene_biotype chr_str sequence";
	file.close();
	int orf_index = 0;
    print_genes(gene_models_shell, output_dir + "candidate_orfs",0, orf_index);
    print_genes(gene_models_shell, output_dir + "candidate_orfs",1, orf_index);
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    cout << "\nprint_genes took: " << elapsed.count() << " s";
    
    cout << "\nprinted files";
}

	else if(runMode=="GenerateTranslationProfile")
	{
		//to run after running STAR 

		string genome_path = getArg(parseArgs, "Genome", "", true);
		string riboseq_files_path = getArg(parseArgs, "Riboseq", "", true);
		string candidate_orfs_path = getArg(parseArgs, "CandidateORFs", "", true);
		string output_dir = getArg(parseArgs, "Output", "", false);
		
		
		int min_length = stoi(getArg(parseArgs, "Min_Length", "25", false));
		int max_length = stoi(getArg(parseArgs, "Max_Length", "35", false));
		int threads = stoi(getArg(parseArgs, "Threads", "1", false));
		float p_site_factor = stof(getArg(parseArgs, "P_Site_Factor", "1.75", false));
		float p_site_distance = stoi(getArg(parseArgs, "P_Site_Distance", "40", false));
		int cutoff = stoi(getArg(parseArgs, "QC_Count", "10000", false)); 
		float required_frame_difference = stof(getArg(parseArgs, "QC_Periodicity", "2.0", false)); 
		bool p_site_old = false;
		
		bool qc_positions = stringToBool(getArg(parseArgs, "QC_Positions", "false", false));
		

		vector<Study> riboseq_studies;
		read_riboseq_studies(riboseq_studies, riboseq_files_path);
		
		int align_start = 0;
		int align_end = riboseq_studies.size()-1;


		if(output_dir.size()!=0){
			create_dir(output_dir);
			
			output_dir+="/";
		}

		
		
		cout<<"\nfile indices: "<<align_start<<" "<<align_end<<" "<<argc<<" "<<riboseq_studies.size();

		map<string, int> chr_labels;
		read_genome_assign_reads(chr_labels, genome_path);



		vector<GeneModel> canonical_orfs;


		//First we read annotated only, read all later on.
		read_genes(canonical_orfs, candidate_orfs_path, true);
		//read_genes(orfs, "orfs");
		//read_genes_yeast(all_orfs, "orfs_comp");
		//read_genes_yeast(orfs, "orfs_comp");

		//read_annotated_genes(canonical_orfs, "orf_coding_all.fasta");
		
		cout << "\norfs read: " << canonical_orfs.size();


		expand_gene_models(canonical_orfs, threads);

		vector<map<int,int>> all_passed_reads_f(chr_labels.size());
		vector<map<int,int>> all_passed_reads_r(chr_labels.size());
		
		
		//Initialize global reads_map variable
		for (int i = align_start; i <= align_end; i++)
		{
		

			cout << "\nprocessing study: " << riboseq_studies[i].path;

			auto time3 = high_resolution_clock::now();

			vector<string> sam_lines;
			read_sam_file(riboseq_studies[i].path,  chr_labels, threads, sam_lines);
			vector<Read> reads(sam_lines.size());

			process_sam_lines(reads, riboseq_studies[i].path, chr_labels, sam_lines, threads);


			cout << "\nreads read: " << reads.size();
			
			auto time4 = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(time4 - time3);
			cout << "Read sam: " << duration.count() << " microseconds" << endl;



			vector<vector<map<int, int>>> reads_map_f;//readlength x chromosome x position
			vector<vector<map<int, int>>> reads_map_r;//readlength x chromosome x position


			reads_map_f = vector<vector<map<int, int>>>(max_length,vector<map<int,int>>(chr_labels.size()));
			reads_map_r = vector<vector<map<int, int>>>(max_length,vector<map<int,int>>(chr_labels.size()));

			map_reads_to_genome(reads_map_f, reads_map_r, reads, chr_labels.size(), threads);


			cout << "\n chr size: " << chr_labels.size();

			
			
			auto time5 = high_resolution_clock::now();
			duration = duration_cast<microseconds>(time5 - time4);
			cout << "Map reads to genome: " << duration.count() << " microseconds" << endl;
		
			reads.clear();
			




			int stop = max_length;
			int start = min_length;

			vector<vector<int>> ann_frames_f(stop-start+1, vector<int>(3,0)); 
			vector<vector<int>> ann_frames_r(stop-start+1, vector<int>(3,0));




			vector<int> p_sites_f(stop-start+1);
			vector<int> p_sites_r(stop-start+1);

			#pragma omp parallel for num_threads(threads)
			for(int j=start; j<=stop; j++){
				int pf = find_p_site(canonical_orfs, reads_map_f[j]);
				int pr = find_p_site_rc(canonical_orfs, reads_map_r[j]);
				
				#pragma omp critical
				{
					p_sites_f[j-start] = pf;
					p_sites_r[j-start] = pr;
				}

			}
			
			for (int i = 0; i < p_sites_f.size(); i++) {
				std::cout << "p_sites_f[" << start + i << "] = " << p_sites_f[i] << std::endl;
				std::cout << "p_sites_r[" << start + i << "] = " << p_sites_r[i] << std::endl;
			}
			
			auto time8 = high_resolution_clock::now();

				
			for (int j = start; j <= stop; j++)
			{

				int p_site_f = p_sites_f[j-start];
				int p_site_r = p_sites_r[j-start]; 
				

				int index = j - start;

				if(p_site_f==-1000 || p_site_r==-1000 || p_site_f-2 != -p_site_r){ //REMEMBER THIS
					cout << "\nRead length: " << to_string(j) << " skipped due to no p-site";
					continue;
				}
				if (quality_control(canonical_orfs, reads_map_f[j], reads_map_r[j], p_site_f, p_site_r, threads, qc_positions) == false){
					cout << "\nRead length: " << to_string(j) << " skipped due to poor quality.";

					continue;
				}


				//assign_reads_to_orfs(all_orfs, reads_map_f[j],reads_map_r[j], ps[index], rps[index], threads);

				accumulate_reads(threads, chr_labels, reads_map_f[j], reads_map_r[j], all_passed_reads_f, all_passed_reads_r, p_site_f, p_site_r);
			}

			
			auto time9 = high_resolution_clock::now();
			duration = duration_cast<microseconds>(time9 - time8);
			cout << "\nMap reads to orfs: " << duration.count() << " microseconds" << endl;

			


		}
		
		canonical_orfs.clear();
		
		print_all_passed_reads(all_passed_reads_f, output_dir + "all_passed_reads_f",0);
		print_all_passed_reads(all_passed_reads_r, output_dir + "all_passed_reads_r",1);
		
		
		
		
		//Detecting Translation
		
		vector<GeneModel> all_orfs;
		read_genes(all_orfs, candidate_orfs_path, false);
		expand_gene_models(all_orfs, threads);
		assign_reads_to_orfs(all_orfs, all_passed_reads_f,all_passed_reads_r, threads);
		print_gene_reads_new(all_orfs, output_dir + "orfs_reads"); // _"+to_string(align_start)  + "_" + to_string(align_end)  + "_" + to_string(cutoff)  + "_" + to_string(required_frame_difference));


		print_translation_calls(all_orfs,  output_dir + "translation_calls");
		outputTracks(all_passed_reads_f, all_passed_reads_r, all_orfs, chr_labels, output_dir);

	}

}


