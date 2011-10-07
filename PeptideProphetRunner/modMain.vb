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
    Private Const PEP_PROPHET_ENZYME As String = "tryptic"

    Private mInputFilePath As String
    Private mOutputFolderPath As String
    Private mMaxRuntimeMinutes As Integer = 0       ' Set to a positive value to abort processing if this many minutes elapses

    Private m_PeptideProphet As PeptideProphetLibrary.PeptideProphet

    Private m_PepProphetThread As System.Threading.Thread
    Private m_PepProphetThreadStart As New System.Threading.ThreadStart(AddressOf StartPeptideProphet)
    Private m_PepProphetRetVal As PeptideProphetLibrary.IPeptideProphet.ProcessStatus


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

        Dim dtStartTime As System.DateTime

        If Not CleanupFilePaths() Then
            Return False
        End If


        Try
            'Call peptide prophet using a separate thread so that we can abort it if it runs too long

            m_PepProphetThread = New System.Threading.Thread(m_PepProphetThreadStart)
            m_PepProphetThread.Start()
            dtStartTime = System.DateTime.UtcNow

            ' Wait 2 seconds
            System.Threading.Thread.Sleep(2000)

            Do While m_PeptideProphet Is Nothing AndAlso System.DateTime.UtcNow.Subtract(dtStartTime).TotalSeconds < 20
                ' Wait some more if the peptide prophet object still isn't instantiated
                System.Threading.Thread.Sleep(1000)
            Loop

            Do While (m_PeptideProphet.Status = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_STARTING) OrElse _
                     (m_PeptideProphet.Status = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_RUNNING)

                System.Threading.Thread.Sleep(3000)

                '' RaiseEvent PeptideProphetRunning("Status = " & m_PeptideProphet.Status.ToString, m_PeptideProphet.PercentComplete)

                If mMaxRuntimeMinutes > 0 AndAlso System.DateTime.UtcNow.Subtract(dtStartTime).TotalMinutes >= mMaxRuntimeMinutes Then
                    Console.WriteLine("Peptide prophet has been running for over " & mMaxRuntimeMinutes.ToString & " minutes; aborting")
                    Console.WriteLine()
                    m_PepProphetThread.Abort()
                    Return False
                End If
            Loop

            If m_PepProphetRetVal = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_ERROR Then
                Console.WriteLine("Peptide prophet returned a non-zero error code: " & m_PepProphetRetVal.ToString)
                Console.WriteLine(m_PeptideProphet.ErrMsg)
                Return False
            End If

            Return True

        Catch ex As System.Exception
            Console.WriteLine("Exception while running peptide prophet: " & ex.Message)
            m_PeptideProphet = Nothing
            Return False
        End Try


    End Function

    Private Function SetOptionsUsingCommandLineParameters(ByVal objParseCommandLine As clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

        Dim strValue As String = String.Empty
        Dim strValidParameters() As String = New String() {"I", "O", "T"}
        Dim intValue As Integer

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

                    If .RetrieveValueForParameter("T", strValue) Then
                        If Integer.TryParse(strValue, intValue) Then
                            If intValue > 0 Then mMaxRuntimeMinutes = intValue
                        End If
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

            Console.WriteLine("This program runs Peptide Prophet on the given Sequest Synopsis file")
            Console.WriteLine()
            Console.WriteLine("Program syntax:" & ControlChars.NewLine & IO.Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location))
            Console.WriteLine(" InputFilePath.txt [/O:OutputFolderPath] [/T:MaxRuntimeMinutes]")
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

    Private Sub StartPeptideProphet()

        ' Initialize peptide prophet
        m_PeptideProphet = New PeptideProphetLibrary.PeptideProphet

        Dim StartParams As PeptideProphetLibrary.InitializationParams
        StartParams = New PeptideProphetLibrary.InitializationParams

        StartParams.InputFileName = mInputFilePath
        StartParams.OutputFolderPath = mOutputFolderPath
        StartParams.Enzyme = PEP_PROPHET_ENZYME
        m_PeptideProphet.Setup(StartParams)

        Try
            'Call peptide prophet
            m_PepProphetRetVal = m_PeptideProphet.Start()

        Catch ex As System.Exception
            Console.WriteLine("Error running Peptide Prophet: " & ex.Message & " " & ex.StackTrace)
            m_PepProphetRetVal = PeptideProphetLibrary.IPeptideProphet.ProcessStatus.PP_ERROR
        End Try
    End Sub
  
End Module
