// PeptideProphetLibrary.h

#pragma once

#include <iostream>
#include <vector>
#include "SequestResult.h"
#include "SequestDiscrimScoreCalculator.h"
#include "MixtureModel.h"
#include <fstream>
#include "GlobalVariable.h"
#include <algorithm>
#include <string>
#include <ctype.h>

using namespace std ;

using namespace System;
using namespace System::IO;
using namespace System::Collections;
using namespace System::Runtime::InteropServices;

namespace PeptideProphetLibrary
{
	public __value struct InitializationParams
	{
		System::String* InputFileName;
		System::String* OutputFolderPath;
		System::String* Enzyme;
		int DebugLevel;
		//ILogger Logger; // Include this if using the prism dll logging
	};

	public __gc class AbortException : public System::Exception
	{
		public:
			AbortException(System::String __gc *s) : System::Exception(s) {}
	};

	// Defines minimum required functionality for classes that will process sequest synopsis results files
	public __gc __interface IPeptideProphet
	{
		__value enum ProcessCheckFile
		{
			PP_IOFileExist,
			PP_IOFileNonexistent
		} ;
		__value enum ProcessResults
		{
			PP_SUCCESS = 0, // Operation succeeded
			PP_FAILURE = -1, // Operation failed
			PP_ABORTED = -2  // PeptideProphet aborted
		};
		__value enum ProcessStatus
		{
			PP_STARTING, // PeptideProphet initialization in progress
			PP_RUNNING,  // PeptideProphet running
			PP_COMPLETE, // PeptideProphet successfully completed
			PP_ERROR,    // There was an error somewhere
			PP_ABORTING  // An ABORT command has been received; plugin shutdown in progress
		};

		// (in) – filename and path to sequest synopsis results file
		__property void set_InputFileName(System::String* value);
		__property System::String* get_InputFileName();
		// (in) – filename and path to results file
		__property void set_OutputFilePath(System::String* value);
		__property System::String* get_OutputFilePath();
		// (in) - Enzyme name - right now only tryptic is supported
		__property void set_Enzyme(System::String* value);
		__property System::String* get_Enzyme();

		__property ProcessCheckFile get_FileStatus() ;

		// Allows calling program to determine if processing succeeded
		__property ProcessResults get_Results();
		// Allows calling program to get current status
		__property ProcessStatus get_Status();
		// Error message describing any errors encountered
		__property System::String* get_ErrMsg();
		// Progress indicator, value between 0 and 100 ** May not need this **
		__property float get_PercentComplete();
		// Allows control of debug information verbosity; 0=minimum, 5=maximum verbosity
		__property int	get_DebugLevel();
		__property void set_DebugLevel(int level);

		//__property void set_Logger(ILogger logger);

		// Initializes parameters. Must be called before executing Start()
		void Setup(InitializationParams InitParams);
		// Check whether or not the input file and output file path exist
		IPeptideProphet::ProcessCheckFile FileCheck() ;
		// Starts the peptide prophet file creation process
		IPeptideProphet::ProcessStatus Start();
		// Aborts peptide prophet file creation
		IPeptideProphet::ProcessResults Abort();
	};



	bool SortSequestResultsByScanChargeXCorrDelcn2Peptide(SequestResult &a, SequestResult &b);

	public __gc class PeptideProphet : public PeptideProphetLibrary::IPeptideProphet
	{
		private:
			void LoadSynopsisFile(char *synopsis_file, std::vector<SequestResult> &vectResults, std::vector<DatasetNumMap> &vecDatasetNumMap);
			char* strCopy(char* orig);
			int PValueCalculate(System::String *synopsis_file, System::String *output_file, System::String *output_file_param, System::String *enzyme);

			// Members
			System::String *input_file_name;
			System::String *output_file_path;
			System::String *enzyme;
			System::String *error_msg;
			IPeptideProphet::ProcessCheckFile file_status ;
			IPeptideProphet::ProcessResults results;
			IPeptideProphet::ProcessStatus status;
			int debug_level;
			float percent_complete;
			bool abort;

		public:
			PeptideProphet(void) {}

			// IPeptideProphet Interface definitions:

			// (in) – filename and path to sequest synopsis results file
			__property void set_InputFileName(System::String* value)
			{
				input_file_name = value;
			}
			__property System::String* get_InputFileName()
			{
				return input_file_name;
			}
			// (in) – filename and path to results file
			__property void set_OutputFilePath(System::String* value)
			{
				output_file_path = value;
			}
			__property System::String* get_OutputFilePath()
			{
				return output_file_path;
			}
			// (in) - Enzyme name - right now only tryptic is supported
			__property void set_Enzyme(System::String* value)
			{
				enzyme = value;
			}
			__property System::String* get_Enzyme()
			{
				return enzyme;
			}

			__property ProcessCheckFile get_FileStatus()
			{
				return file_status ;
			}

			// Allows calling program to determine if processing succeeded
			__property ProcessResults get_Results()
			{
				return results;
			}
			// Allows calling program to get current status
			__property ProcessStatus get_Status()
			{
				return status;
			}
			// Error message describing any errors encountered
			__property System::String* get_ErrMsg()
			{
				return error_msg;
			}
			// Progress indicator, value between 0 and 100 ** May not need this **
			__property float get_PercentComplete()
			{
				return percent_complete;
			}
			// Allows control of debug information verbosity; 0=minimum, 5=maximum verbosity
			__property int	get_DebugLevel()
			{
				return debug_level;
			}
			__property void set_DebugLevel(int level)
			{
				debug_level = level;
			}

			//__property void set_Logger(ILogger logger);

			// Initializes parameters. Must be called before executing Start()
			void Setup(InitializationParams InitParams)
			{
				this->DebugLevel = InitParams.DebugLevel;
				this->Enzyme = InitParams.Enzyme;
				this->InputFileName = InitParams.InputFileName;
				this->OutputFilePath = InitParams.OutputFolderPath;
			}

			IPeptideProphet::ProcessCheckFile FileCheck ()
			{
				if (!System::IO::File::Exists(this->InputFileName))
				{
					throw new AbortException(String::Concat("Input file not found: ", this->InputFileName));
				}

				if (!System::IO::Directory::Exists(this->OutputFilePath))
				{
					throw new AbortException(String::Concat("Output folder not found: ", this->OutputFilePath));
				}

				else
				{
					this->file_status = IPeptideProphet::ProcessCheckFile::PP_IOFileExist ;
					return file_status ;
				}
			}


			// Starts the peptide prophet file creation process
			IPeptideProphet::ProcessStatus Start()
			{
				try
				{
					file_status = FileCheck() ;

					status = IPeptideProphet::ProcessStatus::PP_RUNNING ;
					results = IPeptideProphet::ProcessResults::PP_FAILURE ;

					abort = false;

					String *temp1 ;
					String *temp2 ;
					String *temp3 ;
					int at ;

					temp1 = input_file_name ;
					at = temp1->LastIndexOf(S"\\") ;
					temp2 = temp1->Substring(at+1,temp1->Length-4-at-1) ;
					temp3 = String::Concat(temp2, "_PepProphet") ;
					temp3 = String::Concat(temp3, ".txt") ;

					System::String *output_file
						= System::IO::Path::Combine(this->OutputFilePath, temp3);

					System::String *output_file_param
						= System::IO::Path::Combine(this->OutputFilePath, "PeptideProphet_Coefficients.txt" ) ;

					//System::String *output_file
						//= System::IO::Path::Combine(this->OutputFolderPath, new System::String("peptide_prophet_results.txt"));
					this->PValueCalculate(this->InputFileName, output_file, output_file_param, this->Enzyme);

					status = IPeptideProphet::ProcessStatus::PP_COMPLETE;
					results = IPeptideProphet::ProcessResults::PP_SUCCESS;
				}
				catch(AbortException *ex)
				{
					this->error_msg = ex->Message;
					this->file_status = IPeptideProphet::ProcessCheckFile::PP_IOFileNonexistent;
					this->results = IPeptideProphet::ProcessResults::PP_ABORTED ;
					this->status = IPeptideProphet::ProcessStatus::PP_ERROR;
				}
				catch(System::Exception *ex)
				{
					//status = IPeptideProphet::ProcessStatus::PP_ERROR;

					this->error_msg = ex->Message;
					this->results = IPeptideProphet::ProcessResults::PP_FAILURE;
					this->status = IPeptideProphet::ProcessStatus::PP_ERROR;
				}

				return Status;
			}

			// Aborts peptide prophet file creation
			IPeptideProphet::ProcessResults Abort()
			{
				status = IPeptideProphet::ProcessStatus::PP_ABORTING;
				abort = true;
				return IPeptideProphet::ProcessResults::PP_ABORTED;
			}
	};
}
