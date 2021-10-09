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
	public value struct InitializationParams
	{
		System::String^ InputFileName;
		System::String^ OutputFolderPath;
		System::String^ Enzyme;
		int DebugLevel;
		//ILogger Logger; // Include this if using the prism dll logging
	};

	public ref class AbortException : public System::Exception
	{
		public:
			AbortException(System::String ^s) : System::Exception(s) {}
	};

	// Defines minimum required functionality for classes that will process sequest synopsis results files
	public interface class IPeptideProphet
	{
		enum class ProcessCheckFile
		{
			PP_IOFileExist,
			PP_IOFileNonexistent
		} ;
		enum class ProcessResults
		{
			PP_SUCCESS = 0, // Operation succeeded
			PP_FAILURE = -1, // Operation failed
			PP_ABORTED = -2  // PeptideProphet aborted
		};
		enum class ProcessStatus
		{
			PP_STARTING, // PeptideProphet initialization in progress
			PP_RUNNING,  // PeptideProphet running
			PP_COMPLETE, // PeptideProphet successfully completed
			PP_ERROR,    // There was an error somewhere
			PP_ABORTING  // An ABORT command has been received; plugin shutdown in progress
		};

		// (in) – filename and path to sequest synopsis results file
		property System::String^ InputFileName;
		// (in) – filename and path to results file
		property System::String^ OutputFilePath;
		// (in) - Enzyme name - right now only tryptic is supported
		property System::String^ Enzyme;

		property ProcessCheckFile FileStatus { ProcessCheckFile get(); }

		// Allows calling program to determine if processing succeeded
		property ProcessResults Results { ProcessResults get(); }
		// Allows calling program to get current status
		property ProcessStatus Status { ProcessStatus get(); }
		// Error message describing any errors encountered
		property System::String^ ErrMsg { System::String^ get(); }
		// Progress indicator, value between 0 and 100 ** May not need this **
		property float PercentComplete { float get(); }
		// Allows control of debug information verbosity; 0=minimum, 5=maximum verbosity
		property int	DebugLevel;

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

	public ref class PeptideProphet : public PeptideProphetLibrary::IPeptideProphet
	{
		private:
			void LoadSynopsisFile(char *synopsis_file, std::vector<SequestResult> &vectResults, std::vector<DatasetNumMap> &vecDatasetNumMap);
			char* strCopy(char* orig);
			int PValueCalculate(System::String ^synopsis_file, System::String ^output_file, System::String ^output_file_param, System::String ^enzyme);

			// Members
			System::String ^input_file_name;
			System::String ^output_file_path;
			System::String ^enzyme;
			System::String ^error_msg;
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
			property System::String^ InputFileName
			{
				virtual void set(System::String^ value)
				{
					input_file_name = value;
				}
				virtual System::String^ get()
				{
					return input_file_name;
				}
			}
			// (in) – filename and path to results file
			property System::String^ OutputFilePath
			{
				virtual void set(System::String^ value)
				{
					output_file_path = value;
				}
				virtual System::String^ get()
				{
					return output_file_path;
				}
			}
			// (in) - Enzyme name - right now only tryptic is supported
			property System::String^ Enzyme
			{
				virtual void set(System::String^ value)
				{
					enzyme = value;
				}
				virtual System::String^ get()
				{
					return enzyme;
				}
			}

			property IPeptideProphet::ProcessCheckFile FileStatus
			{
				virtual IPeptideProphet::ProcessCheckFile get()
				{
					return file_status;
				}
			}

			// Allows calling program to determine if processing succeeded
			property IPeptideProphet::ProcessResults Results
			{
				virtual IPeptideProphet::ProcessResults get()
				{
					return results;
				}
			}
			// Allows calling program to get current status
			property IPeptideProphet::ProcessStatus Status
			{
				virtual IPeptideProphet::ProcessStatus get()
				{
					return status;
				}
			}
			// Error message describing any errors encountered
			property System::String^ ErrMsg
			{
				virtual System::String^ get()
				{
					return error_msg;
				}
			}
			// Progress indicator, value between 0 and 100 ** May not need this **
			property float PercentComplete
			{
				virtual float get()
				{
					return percent_complete;
				}
			}
			// Allows control of debug information verbosity; 0=minimum, 5=maximum verbosity
			property int	DebugLevel
			{
				virtual int get()
				{
					return debug_level;
				}
				virtual void set(int level)
				{
					debug_level = level;
				}
			}

			//__property void set_Logger(ILogger logger);

			// Initializes parameters. Must be called before executing Start()
			virtual void Setup(InitializationParams InitParams)
			{
				this->DebugLevel = InitParams.DebugLevel;
				this->Enzyme = InitParams.Enzyme;
				this->InputFileName = InitParams.InputFileName;
				this->OutputFilePath = InitParams.OutputFolderPath;
			}

			virtual IPeptideProphet::ProcessCheckFile FileCheck ()
			{
				if (!System::IO::File::Exists(this->InputFileName))
				{
					throw gcnew AbortException(String::Concat("Input file not found: ", this->InputFileName));
				}

				if (!System::IO::Directory::Exists(this->OutputFilePath))
				{
					throw gcnew AbortException(String::Concat("Output folder not found: ", this->OutputFilePath));
				}

				else
				{
					this->file_status = IPeptideProphet::ProcessCheckFile::PP_IOFileExist ;
					return file_status ;
				}
			}


			// Starts the peptide prophet file creation process
			virtual IPeptideProphet::ProcessStatus Start()
			{
				try
				{
					file_status = FileCheck() ;

					status = IPeptideProphet::ProcessStatus::PP_RUNNING ;
					results = IPeptideProphet::ProcessResults::PP_FAILURE ;

					abort = false;

					String ^temp1 ;
					String ^temp2 ;
					String ^temp3 ;
					int at ;

					temp1 = input_file_name ;
					at = temp1->LastIndexOf("\\") ;
					temp2 = temp1->Substring(at+1,temp1->Length-4-at-1) ;
					temp3 = String::Concat(temp2, "_PepProphet") ;
					temp3 = String::Concat(temp3, ".txt") ;

					System::String ^output_file
						= System::IO::Path::Combine(this->OutputFilePath, temp3);

					System::String ^output_file_param
						= System::IO::Path::Combine(this->OutputFilePath, "PeptideProphet_Coefficients.txt" ) ;

					//System::String ^output_file
						//= System::IO::Path::Combine(this->OutputFolderPath, gcnew System::String("peptide_prophet_results.txt"));
					this->PValueCalculate(this->InputFileName, output_file, output_file_param, this->Enzyme);

					status = IPeptideProphet::ProcessStatus::PP_COMPLETE;
					results = IPeptideProphet::ProcessResults::PP_SUCCESS;
				}
				catch(AbortException ^ex)
				{
					this->error_msg = ex->Message;
					this->file_status = IPeptideProphet::ProcessCheckFile::PP_IOFileNonexistent;
					this->results = IPeptideProphet::ProcessResults::PP_ABORTED ;
					this->status = IPeptideProphet::ProcessStatus::PP_ERROR;
				}
				catch(System::Exception ^ex)
				{
					//status = IPeptideProphet::ProcessStatus::PP_ERROR;

					this->error_msg = ex->Message;
					this->results = IPeptideProphet::ProcessResults::PP_FAILURE;
					this->status = IPeptideProphet::ProcessStatus::PP_ERROR;
				}

				return Status;
			}

			// Aborts peptide prophet file creation
			virtual IPeptideProphet::ProcessResults Abort()
			{
				status = IPeptideProphet::ProcessStatus::PP_ABORTING;
				abort = true;
				return IPeptideProphet::ProcessResults::PP_ABORTED;
			}
	};
}
