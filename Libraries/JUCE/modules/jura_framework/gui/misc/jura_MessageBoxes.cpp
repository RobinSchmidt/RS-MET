#include "rojue_MessageBoxes.h"

void rojue::showMemoryAllocationErrorBox(String sourceFunctionName)
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Memory Allocation Error")), 
    String(T("Memory could not be allocated in: ")) + sourceFunctionName, String(T("OK")) );
}

void rojue::showAudioFileInvalidErrorBox(String fileName)
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Audiofile invalid")), 
    String(T("The file: ")) + fileName 
    + String(T(" is not a supported audiofile.")), String(T("OK")) );
}

void rojue::showAudioFileTooLargeErrorBox()
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Audiofile too large")), 
    String(T("Audiofiles with more than ")) + String(INT_MAX) 
    + String(T(" samples are not supported")), String(T("OK")) );
}

void rojue::showAudioFileWriteErrorBox(String fileName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Audiofile write error")), 
    String(T("The audio-file: ")) + fileName 
    + String(T(" could not be written")), String(T("OK")) );
}

void rojue::showDirectoryCouldNotBeCreatedBox(File directoryThatCouldNotBeCreated)
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Directory could not be created")), 
    String(T("The directory: ")) +  directoryThatCouldNotBeCreated.getFullPathName()
    + String(T(" could not be created")), String(T("OK")) );
}

bool rojue::showDirectoryDeleteWarningBox(String directoryName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Delete directory?")), 
    String(T("The directory: ")) + directoryName 
    + String(T(" will be deleted from the project directory and all clips using a sample from this directory will be removed - proceed?")), 
    String(T("OK")), String(T("Cancel")) );
}

void rojue::showEnterNameErrorBox()
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Please enter a name")), 
    String(T("Please enter a name for the new project in the textfield.")), String(T("OK")) );
}

void rojue::showFileCouldNotBeCopiedBox(const String& fileName, const String& attemptedDirectoryName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("File could not be copied")), 
    String(T("The file: ")) + fileName + String(T(" could not be copied into the directory")) 
    + attemptedDirectoryName, String(T("OK")) );
}

void rojue::showFileCouldNotBeDeletedBox(const String& fileName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("File could not be deleted")), 
    String(T("The file: ")) + fileName + String(T(" could not be deleted")), String(T("OK")) );
}

void rojue::showFileCouldNotBeMovedBox(const String& fileName, const String& attemptedDirectoryName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("File could not be moved")), 
    String(T("The file: ")) + fileName + String(T(" could not be moved into the directory")) 
    + attemptedDirectoryName, String(T("OK")) );
}

void rojue::showFileIsNoDirectoryErrorBox(const File& file)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Not a directory")), 
    String(T("The chosen file: ")) + file.getFileName() 
    + String(T(" is not a directory")), String(T("OK")) );
}

void rojue::showFileNotFoundOrInvalidAudioFileBox(const String& fileName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("File not found or invalid")), 
    String(T("The file: ")) + fileName
    + String(T(" was not found or is not a valid audiofile - clips using this file will be omitted")), 
    String(T("OK")) );
}


void rojue::showInvalidPathNameErrorBox(const String& theInvalidName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Invalid path name")), 
    String(T("The chosen project directory name: ")) + theInvalidName 
    + String(T(" is not a valid name for a directory. Please enter another name.")), String(T("OK")) );
}

bool rojue::showOverwriteAudioFileWarningBox(String fileName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Overwrite audio file?")), 
    String(T("The audio-file: ")) + fileName 
    + String(T(" will be overwritten - proceed?")), String(T("OK")), String(T("Cancel")) );
}

bool rojue::showOverwriteSongFileWarningBox(String fileName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Overwrite song file?")), 
    String(T("The song-file: ")) + fileName 
    + String(T(" will be overwritten - proceed?")), String(T("OK")), String(T("Cancel")) );
}

void rojue::showPleaseSelectProjectFileBox()
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Please select a project file")), 
    String(T("Please select a valid Mixsonic project file to load.")), String(T("OK")) );
}

void rojue::showProjectCouldNotBeCreatedBox(const String& projectName)
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Project could not be created")), 
    String(T("The project: ")) +projectName + String(T(" could not be created")), String(T("OK")));
}

void rojue::showProjectDirectoryAlreadyExistsErrorBox(File existingProjectDirectory)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Directory already exists")), 
    String(T("The chosen project directory: ")) + existingProjectDirectory.getFullPathName() 
    + String(T(" already exists. Please enter another name.")), String(T("OK")) );
}

void rojue::showProjectsParentDirectoryInvalidBox(const String& pathString)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Invalid path for projects")), 
    String(T("The parent directory for the projects: ")) + pathString 
    + String(T(" does not exist or is invalid.")), String(T("OK")) );
}

bool rojue::showSampleDeleteWarningBox(String fileName)
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Delete sample?")), 
    String(T("The sample file: ")) + fileName 
    + String(T(" will be deleted from the project directory and all clips using that sample will be removed - proceed?")), 
    String(T("OK")), String(T("Cancel")) );
}

void rojue::showSampleContentPathInvalidBox(const String& pathString)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Invalid sample content path")), 
    String(T("The path: ")) + pathString 
    + String(T(" for the sample content does not exist or is invalid.")), String(T("OK")) );
}

void rojue::showSettingsFileIsMissingBox()
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Settings file missing")), 
    String(T("The global application settings file: MixsonicSettings.xml was not found in the application directory - Mixsonic will use fallback settings, which are likely to be meaningless. Please re-install to resolve this issue.")), 
    String(T("OK")) );
}

void rojue::showSettingsFileIsInvalidBox()
{
  AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("Settings file missing")), 
    String(T("The global application settings file: MixsonicSettings.xml is corrupted - Mixsonic will use fallback settings, which are likely to be meaningless. Please re-install to resolve this issue.")), 
    String(T("OK")) );
}

bool rojue::showSongUnsavedWarningBox()
{
  return AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String(T("Song not saved!")), 
    String(T("The current song is not saved and will be discarded - proceed?")), 
    String(T("OK")), String(T("Cancel")) );
}

void rojue::showTargetFileAlreadyExistsBox(const String& fileName, const String& attemptedDirectoryName)
{
  return AlertWindow::showMessageBox(AlertWindow::WarningIcon, String(T("File could not be moved")), 
    String(T("The file: ")) + fileName + String(T(" could not be moved into the directory")) 
    + attemptedDirectoryName 
    + String(T(" because there already exists a file with the same name.")), String(T("OK")) );
}
