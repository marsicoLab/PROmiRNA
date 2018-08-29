// ============================================================================
// Functions and routines to read lines from a wib file
// ============================================================================
// Author: Brian Caffrey (caffrey@molgen.mpg.de)

#include<iostream>
#include<fstream>
#include <streambuf>
#include<vector>
#include<map>
#include<string.h>

#define MAX_WIG_VALUE       127
#define BIN_TO_VALUE(bin,lowerLimit,range) \
	lowerLimit+(range*((double)bin/(double)MAX_WIG_VALUE))
//typedef std::map<std::pair<CharString, char>, IntTree<GenomicLocation> > geneStartTrees;

struct WigEntry{
	unsigned int chromStart =-1;
	unsigned int count = -1;
	unsigned int offset = -1;
	double dataRange = 0.0;
	double lowerLimit = 0.0;
	std::string chrom;

};

struct loc{
	unsigned start;
	unsigned end;
	double cons=0.0;
	unsigned count=0;
//	std::vector<unsigned> bedIndices;
	unsigned bedIndex;
};



void readWibToString(const std::string& path, std::string & contents);
void parseWigEntry(const std::vector<std::string>& tokens, const std::string& wibFile, double& cons, int& wigCount, const int& start, const int& end);

/*!
 * @fn Tokenize()
 * @brief breaks a string into a vector of strings via a delimeter
 *
 *
 * @signature	 void Tokenize(str, tokens, delimiters)
 *
 * @param[in]		str			string to be split
 * @param[out]		tokens		vector which substrings will be placed into
 * @param[in]		delimiters	one or more delimiters to be split on
 *
 *
 *
 */

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters)
{
	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

	// Find first "non-delimiter".
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);
	while (std::string::npos != pos || std::string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}


/*!
 * @fn parseWigEntry()
 * @brief Takes a wigEntry and string storing a wib file and calculates the conservation at
 *        the start->end genomic position
 *
 * @signature	 void parseWigEntry(wigEntry, wibFile, cons, wigCount, start, end)
 *
 * @param[in]		wigEntry	The wig location to be found and extracted
 * @param[in]		wibFile		The wib string to be used for conservation calculations
 * @param[in, out]	cons		The conservation values (accumulating)
 * @param[in, out]	wigCount	Counts the number of wig sites which were overlapping
 * @param[in]		start		The start position to count from in the wig
 * @param[in]		end			The end position to count from in the wig
 *
 * Function called by checkWigOverlap.
 *
 */


void parseWigEntry(const std::vector<std::string>& wigEntry, const std::string& wibFile, double& cons, unsigned& wigCount, const int start, const int end){

	unsigned chromStart = std::atoi(wigEntry[2].c_str());
	unsigned  chromEnd = std::atoi(wigEntry[3].c_str());
	unsigned count = std::atoi(wigEntry[6].c_str());
	unsigned  offset = std::atoi(wigEntry[7].c_str());
	double dataRange = std::atof(wigEntry[10].c_str());
	double lowerLimit = std::atof(wigEntry[9].c_str());

	std::string tem1 = wibFile.substr(offset,count);

	unsigned startPos=0;
	if(chromStart<start){
		startPos = start-chromStart;
	}
	unsigned endPos=0;
	if (chromEnd > end)
    {
		endPos = chromEnd-end-1;
	}

	for (unsigned i = startPos; i < (count-endPos); ++i)
	{
			cons += BIN_TO_VALUE(tem1[i], lowerLimit, dataRange);
			wigCount++;
	}
}

/*!
 * @fn readWibToString
 * @brief Reads a wib file, stores it in a std::string
 *
 * @signature	 void readWibToString(path, contents)
 *
 * @param[in]		path		path to the wib file
 * @param[in, out]	contents	the string withing which the wib is stored
 *
 * Function called by getConsBedEntries.
 *
 */

void readWibToString(const std::string& path, std::string & contents){

	std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
	if (in)
	{
		in.seekg(0, std::ios::end);
		contents.resize(in.tellg());
		in.seekg(0, std::ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
	}else{
		throw(errno);
	}

}

/*!
 * @fn isOverlapping
 * @brief Returns true or false if two linear regions overlap
 *
 * @signature	 bool isOverlapping(start1, stop1, start2, stop2)
 *
 * @param[in]	start1	start point of first region
 * @param[in]	stop1	stop point of first region
 * @param[in]	start2	start point of second region
 * @param[in]	stop2	stop point of second region
 *
 * Function called by checkWigOverlap().
 *
 */
bool isOverlapping(const int& start1, const int& stop1, const int& start2, const int& stop2){
	return std::max(start1, start2) < std::min(stop1, stop2);
}


/*!
 * @fn checkWigOverlap
 * @brief Takes the map which stores the bedEntry, searches it to see
 *        if a passed wig entry (toks) overlap and then calls parseWigEntry
 *
 * @signature	 void checkWigOverlap(BedMap, wigEntry, wibFile)
 *
 * @param[in, out]	BedMap			Map(chrom->vector of locations), stores regions to be looked up
 * @param[in]		wigPath			tokens
 * @param[in]		wibPath			std::string  of the wib file to be parsed
 *
 * Function called by readWig().
 *
 * Searches BedMap for regions which overlap with the wigEntry
 *
 */
void checkWigOverlap(std::map<std::string, std::vector<loc> >& BedMap, const std::vector<std::string>& wigEntry, const std::string& wibFile){

	//THIS overlap can be performed much more quickly
	//If the Vector is in order then it should be that it can be done in order
	//perhaps have qaurtiles, loop on from the index of the nearest quartile.

	unsigned wigStart = std::stoi(wigEntry[2].c_str());
	unsigned wigEnd = std::stoi(wigEntry[3].c_str());

	auto chr = BedMap.find(wigEntry[1].c_str());

	if(chr!=BedMap.end()){

		for(unsigned i=0;i<chr->second.size();++i){
				if(isOverlapping(chr->second[i].start, chr->second[i].end, wigStart, wigEnd )){
						parseWigEntry(wigEntry, wibFile, chr->second[i].cons, chr->second[i].count, chr->second[i].start, chr->second[i].end);
				}
		}
	}
}


/*!
 * @fn readWig
 * @brief Reads each line of the wig, calls checkWigOverlap to extract conservations
 *
 *
 * @signature	 void readWig(fileName, wibFile, BedMap)

 *
 * @param[in]		stringcmds		std::string of genomic locations seperated by "\n" then "\t"
 * @param[in]		wigPath			Path to the wig file to be parsed
 * @param[in]		wibPath			Path of the wib file to be parsed
 * @param[in, out]	BedMap			Map(chrom->vector of locations), stores regions to be looked up
 *
 * Function called by getConsBedEntries().
 *
 * Reads a wig line by line, calles chceckWigOverlap
 *
 */
void readWig(const std::string& fileName, const std::string& wibFile, std::map<std::string, std::vector<loc> >& BedMap){

	std::fstream fil1;
	fil1.open(fileName.c_str());
	std::string str1;

	while(getline(fil1, str1)){
		std::vector<std::string> toks;
		Tokenize(str1, toks, "\t");
		checkWigOverlap(BedMap, toks, wibFile);
	}
	fil1.close();
}

/*!
 * @fn getConsBedEntries
 * @brief Takes a string of genomic regions (tab delimited) and overlaps with wig regions
 *        Followed by the extraction of conservation values from wib file which is read
 *        to a string within this function. readWig is called to read the wig entries,
 *        overlap them with the bed regions and then extract the conservation values from
 *        the wib string
 *
 *
 * @signature void getConsBedEntries(stringcmds, wigPath, wibPath, BedMap)
 *
 *
 * @param[in]   stringcmds		std::string of genomic locations seperated by "\n" then "\t"
 * @param[in]   wigPath			Path to the wig file to be parsed
 * @param[in]	wibPath			Path of the wib file to be parsed
 * @param[in, out]	BedMap			Map(chrom->vector of locations), stores regions to be looked up
 *
 * Function called by getAverageConservation()
 *
 *
 */

void getConsBedEntries(const std::string& stringcmds, const std::string& wigPath, const std::string& wibPath, std::map<std::string, std::vector<loc> >& BedMap){

	std::vector<std::string> beds;
	Tokenize(stringcmds, beds, "\n");

	//A map of extreme values from the input beds, quick filter
	std::map<std::string, std::pair<unsigned, unsigned> > chroms;

	//Check they don't already exist, i.e. no repeats, looking up twice
	//	std::map<std::string, unsigned> fullBedMap;

	//Here we could have a quick map construction, removing repeat bed entries
	//The beds should be sorted to make it much quicker but need to check if it's possible

	for(int i=0;i<beds.size();++i){

		std::vector<std::string> fields;
//		auto check = fullBedMap.find(beds[i]);

//		if(check == fullBedMap.end()){
			Tokenize(beds[i], fields, "\t");

			loc tempLoc;
			tempLoc.start = std::stoi(fields[1].c_str());
			tempLoc.end = std::stoi(fields[2].c_str());
			tempLoc.cons = 0.0;
			tempLoc.bedIndex = i;

			BedMap[fields[0]].push_back(tempLoc);
		//	fullBedMap[beds[i]]=1;

		//This section is also commented out on purpose but would be possible in
		//the case of repeated regions which PROmirna shouldn't have

		/*}else{


			unsigned temStart = std::stoi(fields[1].c_str());
			unsigned temEnd= std::stoi(fields[2].c_str());

			auto in = BedMap.find(fields[0]);
			for(int j=0;j<in.size();++j){
				if(in->second[j].start==temStart && in->second[j].end==temEnd){
						in->second[j].push_back(i);
				}
			}

		}*/
	}

	std::string contents;
	readWibToString(wibPath, contents);
	std::vector< std::vector<std::string> > WigVec;
	readWig(wigPath,  contents, BedMap);
}

/*!
 * @fn getAverageConservation
 * @brief Takes a vector of Dataset records and assigns conservation scores (Phastcons)
 *
 * @signature void getAverageConservation(avgConservationScores, dataset, wigPath, wibPath)

 *
 * @param[out]  avgConservationScores	Dataset to store the assigned values.
 * @param[in]   dataset					The variable storing the genomic loci
 * @param[in]	wigPath					Path of the wig file to be parsed
 * @param[in]	wibPath					Path of the wib file to be parsed
 *
 * Function called by create_dataset() to assign conservation score to the CAGE regions
 *
 *
 */
void getAverageConservation(std::vector<double> & avgConservationScores, const std::vector<DatasetRecord> & dataset, const CharString & wigPath, const CharString & wibPath)
{
    std::string stringCmd = "";

    for (auto rec : dataset)
    {
        stringCmd = stringCmd + toCString(rec.location.chr) + "\t" + std::to_string(rec.location.start + 1) + "\t" + std::to_string(rec.location.end) + "\n";
    }


	std::map<std::string, std::vector<loc> > BedMap;
	getConsBedEntries(stringCmd, std::string(toCString(wigPath)), std::string(toCString(wibPath)), BedMap);

	for (auto bid : BedMap)
    {
		for (unsigned i = 0; i < bid.second.size(); ++i)
        {
            if (bid.second[i].count > 0)
                avgConservationScores[bid.second[i].bedIndex] = bid.second[i].cons / bid.second[i].count;
            else
                avgConservationScores[bid.second[i].bedIndex] = 0;
		}
	}
}
