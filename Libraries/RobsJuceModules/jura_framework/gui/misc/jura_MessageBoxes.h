#ifndef jura_MessageBoxes_h
#define jura_MessageBoxes_h

// wrap them into a class rsShowMessageBox such that the boxes can be invoke via the syntax:
// rsShowMessageBox::fileNotFound(fileName), etc.

void JUCE_CALLTYPE showWarningBox(const String& title, const String& message);

/** Opens a warning box reporting a memory allocation error. The sourceFunctionName should be the
name of the function from which the box was summoned. */
void showMemoryAllocationErrorBox(juce::String sourceFunctionName);

/** Opens a warning box reporting an invalid audio file - the passed String should be the name of
the file. */
void showAudioFileInvalidErrorBox(juce::String fileName);

/** Opens a warning box reporting a too large audio file */
void showAudioFileTooLargeErrorBox();

/** Opens a warning box reporting that something went wrong when trying to write an audiofile */
void showAudioFileWriteErrorBox(juce::String fileName);

/** Opens an error box indicating that some directory could not be created. */
void showDirectoryCouldNotBeCreatedBox(juce::File directoryThatCouldNotBeCreated);

/** Opens a warning box reporting that the directory will be removed from the SamplePool and
deleted from the HD, if the action is continued. */
bool showDirectoryDeleteWarningBox(juce::String directoryName);

/** Error-box, indicating that the user has not yet eneterd a name (for a new project). */
void showEnterNameErrorBox();

/** Error-box, indicating that a file could not be copied into an attempted target directory. */
void showFileCouldNotBeCopiedBox(const juce::String& fileName,
  const juce::String& attemptedDirectoryName);

/** Error-box, indicating that a file could not be deleted. */
void showFileCouldNotBeDeletedBox(const juce::String& fileName);

/** Error-box, indicating that a file could not be moved into an attempted target directory. */
void showFileCouldNotBeMovedBox(const juce::String& fileName,
  const juce::String& attemptedDirectoryName);

/** Error-box, indicating that a chosen file should have been a directory but wasn't. */
void showFileIsNoDirectoryErrorBox(const juce::File& file);

/** Error box that the file is either nonexistent or not a valid audiofile. */
void showFileNotFoundOrInvalidAudioFileBox(const juce::String& fileName);

/** Shows an error box, indicating that a given String does not represent a valid path. */
void showInvalidPathNameErrorBox(const juce::String& theInvalidName);

/** Opens a warning box reporting that the audiofile will be overwritten, if the action is
continued. */
bool showOverwriteAudioFileWarningBox(juce::String fileName);

/** Error-box, indicating that the user has not selected a valid project file to load. */
void showPleaseSelectProjectFileBox();

/** Opens an error box indicating that some directory could not be created. */
void showProjectCouldNotBeCreatedBox(const juce::String& projectName);

/** Shows an error box that the project directory that the user wants to create already
exists. */
void showProjectDirectoryAlreadyExistsErrorBox(juce::File existingProjectDirectory);

/** Shows an error box that the project parent directory is invalid. \todo: May be later refined
to a dialog which lets the user choose a new directory. */
void showProjectsParentDirectoryInvalidBox(const juce::String& pathString);

/** Opens a warning box reporting that the songfile will be overwritten, if the action is
continued. */
bool showOverwriteSongFileWarningBox(juce::String fileName);

/** Opens a warning box reporting that the audiofile will be removed from the SamplePool and
deleted from the HD, if the action is continued. */
bool showSampleDeleteWarningBox(juce::String fileName);

/** Shows an error box that the sample content directory is invalid. \todo: May be later refined
to a dialog which lets the user choose a new directory. */
void showSampleContentPathInvalidBox(const juce::String& pathString);

/** Error-box, indicating that the global application settings file is missing. */
void showSettingsFileIsMissingBox();

/** Error-box, indicating that the global application settings file is invalid. */
void showSettingsFileIsInvalidBox();

/** Opens a warning box reporting that the current song is not saved and will be discarded, if
the action is continued. */
bool showSongUnsavedWarningBox();

/** Error-box, indicating that a file could not be moved into an attemted target directory
because in the target directory is already a file with this name. */
void showTargetFileAlreadyExistsBox(const juce::String& fileName,
  const juce::String& attemptedDirectoryName);

#endif 