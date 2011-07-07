Option Strict On

' This program calls Peptide Prophet to process the specified synopsis file
'
' Example command line SynFileName.txt
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started April 6, 2009
' Program renamed from TestPeptideProphetConsole to PeptideProphetRunner on July 7, 2011
'
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 

Module modMain

    Public Const PROGRAM_DATE As String = "July 7, 2011"

    Private mInputFilePath As String
    Private mOutputFolderPath As String

    Public Function Main() As Integer
        ' Returns 0 if no error, 1 if an error

        Dim intReturnCode As Integer
        Dim objParseCommandLine As New clsParseCommandLine
        Dim blnProceed As Boolean

        mInputFilePath = String.Empty

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


    Private Function CleanupFilePaths() As Boolean
        ' Returns True if success, False if failure

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


                Return True
            End If
        Catch ex As System.Exception
            System.Console.WriteLine(ControlChars.NewLine & "Error cleaning up the file paths: " & ControlChars.NewLine & ex.Message)
        End Try

        Return False

    End Function

    Private Function GetAppFolderPath() As String
        ' Could use Application.StartupPath, but .GetExecutingAssembly is better
        Return System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location)
    End Function

    Private Function GetAppVersion() As String
        'Return System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")"

        Return System.Reflection.Assembly.GetExecutingAssembly.GetName.Version.ToString & " (" & PROGRAM_DATE & ")"
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

        Dim strValue As String = String.Empty
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
            Console.WriteLine("Error parsing the command line parameters: " & ControlChars.NewLine & ex.Message)
        End Try

        Return False

    End Function

    Private Sub ShowProgramHelp()

        Try

            Console.WriteLine("TThis program runs Peptide Prophet on the given Sequest Synopsis file")
            Console.WriteLine()
            Console.WriteLine("Program syntax:" & ControlChars.NewLine & IO.Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location))
            Console.WriteLine(" InputFilePath.txt [/O:OutputFolderPath]")
            Console.WriteLine()

            Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2003")
            Console.WriteLine("Version: " & GetAppVersion())
            Console.WriteLine()

            Console.WriteLine("E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com")
            Console.WriteLine("Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/")
            Console.WriteLine()

            ' Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
            System.Threading.Thread.Sleep(750)

        Catch ex As System.Exception
            Console.WriteLine("Error displaying the program syntax: " & ex.Message)
        End Try

    End Sub

End Module
