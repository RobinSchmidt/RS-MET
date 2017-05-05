

// maybe move to header file - at least a declaration
void JUCE_CALLTYPE showWarningBox(const String& title, const String& message) // rename to showWarningBox
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, title, message, "OK");
}


void showMemoryAllocationErrorBox(String sourceFunctionName)
{
  showWarningBox("Memory Allocation Error", 
    "Memory could not be allocated in: " + sourceFunctionName);
}

void showAudioFileInvalidErrorBox(String fileName)
{
  showWarningBox("Audiofile invalid", 
    "The file: " + fileName + " is not a supported audiofile.");
}

void showAudioFileTooLargeErrorBox()
{
  showWarningBox("Audiofile too large",  
    "Audiofiles with more than " + String(INT_MAX) + " samples are not supported");
}

void showAudioFileWriteErrorBox(String fileName)
{
  showWarningBox("Audiofile write error", 
    "The audio-file: " + fileName + " could not be written");
}

void showDirectoryCouldNotBeCreatedBox(File directoryThatCouldNotBeCreated)
{
  showWarningBox("Directory could not be created", 
    "The directory: " +  directoryThatCouldNotBeCreated.getFullPathName()
    + " could not be created");
}

void showEnterNameErrorBox()
{
  showWarningBox("Please enter a name", 
    "Please enter a name for the new project in the textfield.");
}

void showFileCouldNotBeCopiedBox(const String& fileName, const String& attemptedDirectoryName)
{
  showWarningBox("File could not be copied", 
    "The file: " + fileName + " could not be copied into the directory") + attemptedDirectoryName);
}

void showFileCouldNotBeDeletedBox(const String& fileName)
{
  showWarningBox("File could not be deleted", 
    "The file: " + fileName + " could not be deleted");
}

void showFileCouldNotBeMovedBox(const String& fileName, const String& attemptedDirectoryName)
{
  showWarningBox("File could not be moved", 
    "The file: " + fileName + " could not be moved into the directory" + attemptedDirectoryName);
}

void showWarningBox(const File& file)
{
  showWarningBox(AlertWindow::WarningIcon, String(T("Not a directory")), 
    String(T("The chosen file: ")) + file.getFileName() 
    + String(T(" is not a directory")), String(T("OK")) );
}

void showWarningBox(const String& fileName)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("File not found or invalid")), 
    String(T("The file: ")) + fileName
    + String(T(" was not found or is not a valid audiofile - clips using this file will be omitted")), 
    String(T("OK")) );
}

void showWarningBox(const String& theInvalidName)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Invalid path name")), 
    String(T("The chosen project directory name: ")) + theInvalidName 
    + String(T(" is not a valid name for a directory. Please enter another name.")), String(T("OK")) );
}


void showPleaseSelectProjectFileBox()
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Please select a project file")), 
    String(T("Please select a valid Mixsonic project file to load.")), String(T("OK")) );
}

void showProjectCouldNotBeCreatedBox(const String& projectName)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Project could not be created")), 
    String(T("The project: ")) +projectName + String(T(" could not be created")), String(T("OK")));
}

void showProjectDirectoryAlreadyExistsErrorBox(File existingProjectDirectory)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Directory already exists")), 
    String(T("The chosen project directory: ")) + existingProjectDirectory.getFullPathName() 
    + String(T(" already exists. Please enter another name.")), String(T("OK")) );
}

void showProjectsParentDirectoryInvalidBox(const String& pathString)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Invalid path for projects")), 
    String(T("The parent directory for the projects: ")) + pathString 
    + String(T(" does not exist or is invalid.")), String(T("OK")) );
}

void showSampleContentPathInvalidBox(const String& pathString)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Invalid sample content path")), 
    String(T("The path: ")) + pathString 
    + String(T(" for the sample content does not exist or is invalid.")), String(T("OK")) );
}

void showSettingsFileIsMissingBox()
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Settings file missing")), 
    String(T("The global application settings file: MixsonicSettings.xml was not found in the application directory - Mixsonic will use fallback settings, which are likely to be meaningless. Please re-install to resolve this issue.")), 
    String(T("OK")) );
}

void showSettingsFileIsInvalidBox()
{
  showMessageBox(AlertWindow::WarningIcon, String(T("Settings file missing")), 
    String(T("The global application settings file: MixsonicSettings.xml is corrupted - Mixsonic will use fallback settings, which are likely to be meaningless. Please re-install to resolve this issue.")), 
    String(T("OK")) );
}

void showTargetFileAlreadyExistsBox(const String& fileName, const String& attemptedDirectoryName)
{
  showMessageBox(AlertWindow::WarningIcon, String(T("File could not be moved")), 
    String(T("The file: ")) + fileName + String(T(" could not be moved into the directory")) 
    + attemptedDirectoryName 
    + String(T(" because there already exists a file with the same name.")), String(T("OK")) );
}


// OK/Cancel boxes:

bool showDirectoryDeleteWarningBox(String directoryName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Delete directory?")), 
    String(T("The directory: ")) + directoryName 
    + String(T(" will be deleted from the project directory and all clips using a sample from this directory will be removed - proceed?")), 
    String(T("OK")), String(T("Cancel")) );
}

bool showOverwriteAudioFileWarningBox(String fileName)
{
  showOkCancelBox(AlertWindow::WarningIcon, String(T("Overwrite audio file?")), 
    String(T("The audio-file: ")) + fileName 
    + String(T(" will be overwritten - proceed?")), String(T("OK")), String(T("Cancel")) );
}

bool showOverwriteSongFileWarningBox(String fileName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Overwrite song file?")), 
    String(T("The song-file: ")) + fileName 
    + String(T(" will be overwritten - proceed?")), String(T("OK")), String(T("Cancel")) );
}

bool showSampleDeleteWarningBox(String fileName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Delete sample?")), 
    String(T("The sample file: ")) + fileName 
    + String(T(" will be deleted from the project directory and all clips using that sample will be removed - proceed?")), 
    String(T("OK")), String(T("Cancel")) );
}

bool showSongUnsavedWarningBox()
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Song not saved!")), 
    String(T("The current song is not saved and will be discarded - proceed?")), 
    String(T("OK")), String(T("Cancel")) );
}

