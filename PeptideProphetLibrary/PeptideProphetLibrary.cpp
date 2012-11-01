// This is the main DLL file.

#include "PeptideProphetLibrary.h"
#include "DatasetNumMap.h"

namespace PeptideProphetLibrary
{

	char* PeptideProphet::strCopy(char* orig) 
	{
		char* output = new char[strlen(orig)+1];
		strcpy(output, orig);
		output[strlen(orig)] = 0;
		return output;
	}

	bool SortSequestResultsByScanChargeXCorrDelcn2Peptide(SequestResult &a, SequestResult &b)
	{
		if (a.dataset_num_ < b.dataset_num_)
			return true ;
		if (a.dataset_num_ > b.dataset_num_)
			return false ;

		if(a.ScanNumber < b.ScanNumber)
			return true ; 
		if (a.ScanNumber > b.ScanNumber)
			return false ; 

		if (a.charge_ < b.charge_)
			return true ; 
		if (a.charge_ > b.charge_)
			return false ;

		if (a.xcorr_ < b.xcorr_)
			return true ; 
		if (a.xcorr_ > b.xcorr_)
		return false ;

		if (a.delta_ < b.delta_)
			return true ; 
		if (a.delta_ > b.delta_)
			return false ;

		if (a.peptide_ < b.peptide_)
			return true ; 
		if (a.peptide_ > b.peptide_)
			return false ;

		return false ; 
	}

	void PeptideProphet::LoadSynopsisFile(char *synopsis_file, std::vector<SequestResult> &vectResults, std::vector<DatasetNumMap> &vecDatasetNumMap)
	{
		//Xiuxia 05/16/2006
		int temp;
		int dataset_num, Scan, NumScans, Charge, RankSp, RankXc, PassFilt, NTT;
		double MH, XCorr, DeltaCn, Sp, MO, DeltaCn2, DelM, XcRatio, MScore;
		char Reference[256];
		char peptide[256] ;

		DatasetNumMap oneMap ;

		dataset_num = 0;
		Scan = 0;
		NumScans = 0;
		Charge = 0;
		MH = 0.0;
		XCorr = 0.0;
		DeltaCn = 0.0;
		Sp = 0.0;
		MO = 0;
		DeltaCn2 = 0.0;
		RankSp = 0;
		RankXc = 0;
		DelM = 0.0;
		XcRatio = 0.0;
		PassFilt = 0;
		MScore = 0.0;
		NTT = 0;

		SequestResult result ; 
		
		if(strlen(synopsis_file)==0)
		{
			System::String *err = new System::String(" Synopsis file path is null");
			throw new System::Exception(err);
		}

		std::ifstream fhtml(synopsis_file);

		if(! fhtml) {
			//std::cerr << " Could not open synopsis file " << synopsis_file << std::endl;
			std::cout << " Could not open synopsis file, file not found" << synopsis_file << "\n";
			System::String *err = new System::String(" Could not open synopsis file, file not found");
			throw new System::Exception(err);
			//exit(1);
		}

		const int MAX_BUFFER_SIZE = 2048 ;
		char data[MAX_BUFFER_SIZE];

		int numLoaded = 0 ; 
		int numSkipped = 0 ;
		float tolerance = 0.00001 ; 
		//Xiuxia, if this tolerance is too small, then equal numbers might be considered unequal, 06/06/2006

		char testLetter ;
		string teststring ;

		fhtml.getline(data, MAX_BUFFER_SIZE) ;	//Xiuxia, read the first line of the input file
		teststring = data ;

		if (teststring.length() != 0) {//Xiuxia, test if the first line is empty

			testLetter = data[0] ;
			if (isdigit(testLetter) == 0) {
				// the first line is header and skip

				;
			}
			else {
				// the first line is data and read it in

				temp = sscanf(data, "%ld %d %d %d %lf %lf %lf %lf %s %lf %s %lf %d %d %lf %lf %d %lf %d",
					&dataset_num, &Scan, &NumScans, &Charge, &MH, &XCorr, &DeltaCn, &Sp, Reference, &MO, peptide, &DeltaCn2, 
					&RankSp, &RankXc, &DelM, &XcRatio, &PassFilt, &MScore, &NTT) ;

				result.dataset_num_ = dataset_num ; 
				result.ScanNumber = Scan;
				result.charge_ = Charge;
				result.xcorr_ = XCorr;
				result.delta_ = DeltaCn2;
				result.rank_ = RankSp;
				result.massdiff_ = DelM;
				result.peptide_ = peptide ;
				result.protein_ = Reference ; 
				result.degen_ = (MO > 0) ; 

				// now add result to vector.
				if (Charge <= numCharge)
				{
					vectResults.push_back(result) ; 
					numLoaded++ ; 
				}
			}
		}

		
		//Xiuxia, read in the rest of the lines
		while(fhtml.getline(data, MAX_BUFFER_SIZE)) 
		{		
			// read line by line and set every member variable of SequestResult (other than discriminant score)
			// for each line read. HERE

			teststring = data ;
			if (teststring.length() != 0) {

				temp = sscanf(data, "%ld %d %d %d %lf %lf %lf %lf %s %lf %s %lf %d %d %lf %lf %d %lf %d",
					&dataset_num, &Scan, &NumScans, &Charge, &MH, &XCorr, &DeltaCn, &Sp, Reference, &MO, peptide, &DeltaCn2, 
					&RankSp, &RankXc, &DelM, &XcRatio, &PassFilt, &MScore, &NTT) ;

				result.dataset_num_ = dataset_num ; 
				result.ScanNumber = Scan;
				result.charge_ = Charge;
				result.xcorr_ = XCorr;
				result.delta_ = DeltaCn2;
				result.rank_ = RankSp;
				result.massdiff_ = DelM;
				result.peptide_ = peptide ;
				result.protein_ = Reference ; 
				result.degen_ = (MO > 0) ; 

				// now add result to vector.
				if (Charge <= numCharge)
				{
					vectResults.push_back(result) ; 
					numLoaded++ ; 
				}
				else
				{
					numSkipped++ ;
				}
			}
		}

		std::cout << " Loaded  " << numLoaded << " results\n" ;
		std::cout << " Skipped " << numSkipped << " results where charge is " << (numCharge+1) << "+ or higher\n" ;
		fhtml.close() ; 
		
		sort(vectResults.begin(), vectResults.end(), SortSequestResultsByScanChargeXCorrDelcn2Peptide) ; 

		if (vectResults.size() == 0)
			return ; 

		std::vector<SequestResult> copyVect ; 
		copyVect.push_back(vectResults[0]) ; 

		int masterDatasetNum ;

		masterDatasetNum = vectResults[0].dataset_num_ ;
		oneMap.dataset_num_start = vectResults[0].dataset_num_ ;
		oneMap.dataset_num_other = vectResults[0].dataset_num_ ;

		vecDatasetNumMap.push_back(oneMap) ;

		for (int index = 1 ; index < vectResults.size() ; index++)
		{ 

			int cDatasetNum = vectResults[index].dataset_num_ ;
			int cScan = vectResults[index].ScanNumber ; 
			int clast_Scan = vectResults[index-1].ScanNumber ; 
			int cCharge = vectResults[index].charge_ ; 
			int clast_Charge = vectResults[index-1].ScanNumber ; 
			double cXCorr = vectResults[index].xcorr_ ; 
			double clast_XCorr = vectResults[index-1].xcorr_ ; 
			double cDeltaCn2 = vectResults[index].delta_ ; 
			double clast_DeltaCn2 = vectResults[index-1].delta_ ; 
			string cPeptide = vectResults[index].peptide_ ;
			string clast_Peptide = vectResults[index-1].peptide_ ;

			if ((cScan != clast_Scan) || (cCharge != clast_Charge) || 
				abs(cXCorr - clast_XCorr) >= tolerance || 
				abs(cDeltaCn2 - clast_DeltaCn2) >= tolerance || 
				(cPeptide != clast_Peptide)) 

			{
				copyVect.push_back(vectResults[index]) ;

				masterDatasetNum = vectResults[index].dataset_num_;

				oneMap.dataset_num_start = masterDatasetNum ;
				oneMap.dataset_num_other = masterDatasetNum ;

				vecDatasetNumMap.push_back(oneMap) ;

			}
			else
			{
				oneMap.dataset_num_start = masterDatasetNum ;
				oneMap.dataset_num_other = cDatasetNum ;

				vecDatasetNumMap.push_back (oneMap) ;
			}

		}

		vectResults.clear() ; 
		vectResults.insert(vectResults.begin(), copyVect.begin(), copyVect.end()) ; 
		copyVect.clear() ; 

		std::cout << " Total number of unique results to be processed = " << vectResults.size() << "\n\n" ;
	}



	int PeptideProphet::PValueCalculate(System::String *synopsis_file, System::String *output_file, System::String *output_file_param, System::String *p_enzyme)
	{	
		std::cout << "\nPeptideProphet v. 1.0 by A.Keller 11.7.02 ISB \n" ;
		std::cout << "Modified by Xiuxia Du and Deep Jaitly, June 21, 2006\n" ;
		std::cout << "Last updated November 1, 2012\n" ;

		//if (argc != 3 && argc != 4)
		//	Usage(argc, argv) ; 

		std::vector<SequestResult> vectResults ; 
		std::vector<DatasetNumMap> vecDatasetNumMap ;

		// load the synopsis file into a vector of results.
		char* c_synopsis_file;
		try
		{
			c_synopsis_file = (char*)(void*)Marshal::StringToHGlobalAnsi(synopsis_file);

			std::cout << "\nLoading synopsis file, " << c_synopsis_file << "\n";
			std::cout.flush();

			LoadSynopsisFile(c_synopsis_file, vectResults, vecDatasetNumMap) ; 

			Marshal::FreeHGlobal(c_synopsis_file);	
		}
		catch(System::Exception *ex)
		{
			Marshal::FreeHGlobal(c_synopsis_file);

			this->error_msg = ex->Message;
			this->results = IPeptideProphet::ProcessResults::PP_FAILURE;
			this->status = IPeptideProphet::ProcessStatus::PP_ERROR;

			std::cout << "Error: " <<  ex->Message << "\n";
			std::cout.flush();

			return 1;
		}

		// now calculate all the discriminant scores. 
		ScoreCalculator* calc = NULL;
		Boolean2 exclude = False;
		Boolean2 windows = False;
		Boolean2 massd = False;
		Boolean2 modify_deltas = False;
		Boolean2 maldi = False;
		int MAX_NUM_ITERS = 500;
		Boolean2 icat = False;
		Boolean2 glyc = False;
		Boolean2 mascot = False;
		Boolean2 qtof = False;

		std::cout << "Processing\n";
		std::cout.flush();

		char *c_enzyme;
		char* c_output_file;
		char *c_output_file_param ;

		try
		{
			if(abort)		
				throw new AbortException("Aborted");

			// Calculate SequestDiscrimScore
			c_enzyme = (char*)(void*)Marshal::StringToHGlobalAnsi(p_enzyme); 			

			// Note: Matt Monroe and Deep Jaitly tried to consolidate the two calls to SequestDiscrimScoreCalculator by
			// updating c_enzyme to be "tryptic" if p_enzyme = "", but were not successful, most likely due to a memory
			// allocation error elsewhere in this application
			if(strlen(c_enzyme) == 0)
				calc = new SequestDiscrimScoreCalculator(vectResults, exclude, windows, massd, modify_deltas, maldi, "tryptic");
			else				
				calc = new SequestDiscrimScoreCalculator(vectResults, exclude, windows, massd, modify_deltas, maldi, c_enzyme);				

			if(abort)
				throw new AbortException("Aborted");

			// See above comment concerning c_enzyme
			MixtureModel* mixmodel ;
			if(strlen(c_enzyme) == 0)
				mixmodel = new MixtureModel(vectResults, MAX_NUM_ITERS, icat, glyc, massd, mascot, qtof, "tryptic");
			else
				mixmodel = new MixtureModel(vectResults, MAX_NUM_ITERS, icat, glyc, massd, mascot, qtof, c_enzyme);
			Marshal::FreeHGlobal(c_enzyme);

			SequestDiscrimFunction* SequestDiscrimFcn = new SequestDiscrimFunction(0) ;

			if(abort)		
				throw new AbortException("Aborted");

			// Write results
			c_output_file = (char*)(void*)Marshal::StringToHGlobalAnsi(output_file);
			c_output_file_param = (char*)(void*)Marshal::StringToHGlobalAnsi(output_file_param);

			//mixmodel->writeResults(c_output_file) ;

			std::cout << "Saving results" << std::endl;

			mixmodel->writeResultsOrdered(c_output_file, vecDatasetNumMap) ;

			SequestDiscrimFcn->writeCoef(c_output_file_param) ; 
			
			Marshal::FreeHGlobal(c_output_file);	
		}
		catch(AbortException *ex)
		{
			Marshal::FreeHGlobal(c_enzyme);
			Marshal::FreeHGlobal(c_output_file);

			this->error_msg = ex->Message;
			this->results = IPeptideProphet::ProcessResults::PP_ABORTED;
			this->status = IPeptideProphet::ProcessStatus::PP_COMPLETE;
			return 1;
		}
		//catch(System::Exception *ex)
		//{
		//	Marshal::FreeHGlobal(c_enzyme);
		//	Marshal::FreeHGlobal(c_output_file);			

		//	this->error_msg = ex->Message;
		//	this->results = IPeptideProphet::ProcessResults::PP_FAILURE;
		//	this->status = IPeptideProphet::ProcessStatus::PP_ERROR;

		//	std::cout << results << "\t" << status << std::endl ;

		//	return 1;
		//}

		return 0 ; 
	}

}