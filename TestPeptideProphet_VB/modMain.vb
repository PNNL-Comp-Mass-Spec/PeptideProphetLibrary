Option Strict On

' This program calls Peptide Prophet to process the specified synopsis file
'
' Example command line /I:SynFileName.txt
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started April 6, 2009
'
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'
' Notice: This computer software was prepared by Battelle Memorial Institute, 
' hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
' Department of Energy (DOE).  All rights in the computer software are reserved 
' by DOE on behalf of the United States Government and the Contractor as 
' provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
' SOFTWARE.  This notice including this sentence must appear on any copies of 
' this computer software.


Module modMain

    Public Const PROGRAM_DATE As String = "July 7, 2011"

    Private mInputFilePath As String
    Private mOutputFolderPath As String
    Private mColumnNumber As Integer
    Private mColumnDelimiter As Char

    Private Function CleanupFilePaths() As Boolean
        ' Returns True if success, False if failure

        Dim intCharLoc As Integer
        Dim strParentFolderPath As String

        Try
            ' Make sure mInputFilePath points to a valid file
            If Not System.IO.File.Exists(mInputFilePath) Then
                System.Console.WriteLine(ControlChars.NewLine & "Error -- Input file not found: " & ControlChars.NewLine & mInputFilePath)
                CleanupFilePaths = False
            Else
                If mOutputFolderPath Is Nothing OrElse mOutputFolderPath.Length = 0 Then
                    ' Define mOutputFolderPath based on mInputFilePath
                    Dim ioFileInfo As System.IO.FileInfo
                    ioFileInfo = New System.IO.FileInfo(mInputFilePath)
                    mOutputFolderPath = ioFileInfo.DirectoryName
                End If


                CleanupFilePaths = True
            End If
        Catch ex As System.Exception
            System.Console.WriteLine(ControlChars.NewLine & "Error cleaning up the file paths: " & ControlChars.NewLine & ex.Message)
        End Try

    End Function

    Private Function GetAppFolderPath() As String
        ' Could use Application.StartupPath, but .GetExecutingAssembly is better
        Return System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location)
    End Function

    Public Function Main() As Integer
        ' Returns 0 if no error, 1 if an error

        Dim intReturnCode As Integer
        Dim objParseCommandLine As New clsParseCommandLine
        Dim blnProceed As Boolean

        mColumnNumber = 1
        mInputFilePath = String.Empty
        mOutputFolderPath = String.Empty
        mColumnDelimiter = ControlChars.Tab

        Try
            blnProceed = False
            If objParseCommandLine.ParseCommandLine Then
                If SetOptionsUsingCommandLineParameters(objParseCommandLine) Then blnProceed = True
            End If

            If Not blnProceed OrElse _
               objParseCommandLine.NeedToShowHelp OrElse _
               objParseCommandLine.ParameterCount + objParseCommandLine.NonSwitchParameterCount = 0 OrElse _
               mInputFilePath.Length = 0 Then
                ShowProgramHelp()
                intReturnCode = 1
            Else
                If RunPeptideProphet() Then
                    intReturnCode = 0
                Else
                    intReturnCode = 1
                End If
            End If

        Catch ex As System.Exception
            System.Console.WriteLine(ControlChars.NewLine & "Error occurred in modMain->Main: " & ControlChars.NewLine & ex.Message)
            intReturnCode = -1
        End Try

        Return intReturnCode

    End Function

    Private Function RunPeptideProphet() As Boolean

        Const PEP_PROPHET_ENZYME As String = "tryptic"

        If Not CleanupFilePaths() Then
            Return False
        End If

        'Call peptide prophet

        Dim p As PeptideProphetLibrary.PeptideProphet = New PeptideProphetLibrary.PeptideProphet
        Dim StartParams As New PeptideProphetLibrary.InitializationParams
        StartParams.InputFileName = mInputFilePath
        StartParams.OutputFolderPath = mOutputFolderPath
        StartParams.Enzyme = PEP_PROPHET_ENZYME

        p.Setup(StartParams)
        Dim RetVal As PeptideProphetLibrary.IPeptideProphet.ProcessStatus

        Try
            RetVal = p.Start()

            Do While (p.Status = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_STARTING) OrElse _
                     (p.Status = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_RUNNING)
                System.Threading.Thread.Sleep(2000)
            Loop

            If RetVal = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_ERROR Then
                Console.WriteLine("Error while running peptide prophet: " & p.ErrMsg)
                Return False
            End If

            Return True

        Catch ex As System.Exception
            Console.WriteLine("Exception while running peptide prophet: " & ex.Message)
            p = Nothing
            Return False
        End Try


    End Function

    Private Function SetOptionsUsingCommandLineParameters(ByVal objParseCommandLine As clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

        Dim strValue As String
        Dim strValidParameters() As String = New String() {"I", "O"}

        Try
            ' Make sure no invalid parameters are present
            If objParseCommandLine.InvalidParametersPresent(strValidParameters) Then
                Return False
            Else
                With objParseCommandLine
                    ' Query objParseCommandLine to see if various parameters are present
                    If .RetrieveValueForParameter("I", strValue) Then
                        mInputFilePath = strValue
                    ElseIf .NonSwitchParameterCount > 0 Then
                            mInputFilePath = .RetrieveNonSwitchParameter(0)
                    End If

                    If .RetrieveValueForParameter("O", strValue) Then
                        mOutputFolderPath = strValue
                    ElseIf .NonSwitchParameterCount >= 2 Then
                        mOutputFolderPath = .RetrieveNonSwitchParameter(1)
                    End If
                End With

                Return True
            End If

        Catch ex As System.Exception
            System.Console.WriteLine(ControlChars.NewLine & "Error parsing the command line parameters: " & ControlChars.NewLine & ex.Message)
        End Try

    End Function

    Private Sub ShowProgramHelp()

        Dim strSyntax As String
        Dim ioPath As System.IO.Path

        Try

            strSyntax = "This program runs Peptide Prophet on the given Synopsis file" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Program syntax:" & ControlChars.NewLine & ioPath.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location)
            strSyntax &= " /I:SynFilePath.txt /O:OutputFolderPath" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2003" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "This is version " & System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com" & ControlChars.NewLine
            strSyntax &= "Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Licensed under the Apache License, Version 2.0; you may not use this file except in compliance with the License.  "
            strSyntax &= "You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Notice: This computer software was prepared by Battelle Memorial Institute, "
            strSyntax &= "hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the "
            strSyntax &= "Department of Energy (DOE).  All rights in the computer software are reserved "
            strSyntax &= "by DOE on behalf of the United States Government and the Contractor as "
            strSyntax &= "provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY "
            strSyntax &= "WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS "
            strSyntax &= "SOFTWARE.  This notice including this sentence must appear on any copies of "
            strSyntax &= "this computer software." & ControlChars.NewLine

            System.Console.WriteLine(ControlChars.NewLine & strSyntax)

        Catch ex As System.Exception
            MsgBox("Error displaying the program syntax: " & ControlChars.NewLine & ex.Message, MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
        End Try

    End Sub

End Module
