/*
  ==============================================================================

   This file is part of the JUCE library.
   Copyright (c) 2017 - ROLI Ltd.

   JUCE is an open source library subject to commercial or open-source
   licensing.

   By using JUCE, you agree to the terms of both the JUCE 5 End-User License
   Agreement and JUCE 5 Privacy Policy (both updated and effective as of the
   27th April 2017).

   End User License Agreement: www.juce.com/juce-5-licence
   Privacy Policy: www.juce.com/juce-5-privacy-policy

   Or: You may also use this code under the terms of the GPL v3 (see
   www.gnu.org/licenses).

   JUCE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY, AND ALL WARRANTIES, WHETHER
   EXPRESSED OR IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR PURPOSE, ARE
   DISCLAIMED.

  ==============================================================================
*/

#include "../Application/jucer_Headers.h"
#include "../ProjectSaving/jucer_ProjectSaver.h"
#include "../ProjectSaving/jucer_ProjectExport_Xcode.h"
#include "../Application/jucer_Application.h"

//==============================================================================
ModuleDescription::ModuleDescription (const File& folder)
   : moduleFolder (folder),
     moduleInfo (parseJUCEHeaderMetadata (getHeader()))
{
}

File ModuleDescription::getHeader() const
{
    if (moduleFolder != File())
    {
        const char* extensions[] = { ".h", ".hpp", ".hxx" };

        for (auto e : extensions)
        {
            auto header = moduleFolder.getChildFile (moduleFolder.getFileName() + e);

            if (header.existsAsFile())
                return header;
        }
    }

    return {};
}

StringArray ModuleDescription::getDependencies() const
{
    auto deps = StringArray::fromTokens (moduleInfo ["dependencies"].toString(), " \t;,", "\"'");
    deps.trim();
    deps.removeEmptyStrings();
    return deps;
}

//==============================================================================
static bool tryToAddModuleFromFolder (const File& path, ModuleIDAndFolderList& list)
{
    ModuleDescription m (path);

    if (m.isValid())
    {
        list.push_back ({ m.getID(), path });
        return true;
    }

    return false;
}

static void addAllModulesInSubfoldersRecursively (const File& path, int depth, ModuleIDAndFolderList& list)
{
    if (depth > 0)
    {
        for (DirectoryIterator iter (path, false, "*", File::findDirectories); iter.next();)
        {
            if (auto* job = ThreadPoolJob::getCurrentThreadPoolJob())
                if (job->shouldExit())
                    return;

            auto childPath = iter.getFile();

            if (! tryToAddModuleFromFolder (childPath, list))
                addAllModulesInSubfoldersRecursively (childPath, depth - 1, list);
        }
    }
}

static void addAllModulesInFolder (const File& path, ModuleIDAndFolderList& list)
{
    if (! tryToAddModuleFromFolder (path, list))
    {
        int subfolders = 3;
        addAllModulesInSubfoldersRecursively (path, subfolders, list);
    }
}

static void sort (ModuleIDAndFolderList& listToSort)
{
    std::sort (listToSort.begin(), listToSort.end(), [] (const ModuleIDAndFolder& m1, const ModuleIDAndFolder& m2)
                                                     {
                                                         return m1.first.compareIgnoreCase (m2.first) < 0;
                                                     });
}

//==============================================================================
struct ModuleScannerJob    : public ThreadPoolJob
{
    ModuleScannerJob (const Array<File>& paths, std::function<void (const ModuleIDAndFolderList&)>&& callback)
        : ThreadPoolJob ("ModuleScannerJob"),
          pathsToScan (paths),
          completionCallback (std::move (callback))
    {
    }

    JobStatus runJob() override
    {
        ModuleIDAndFolderList list;

        for (auto& p : pathsToScan)
            addAllModulesInFolder (p, list);

        if (! shouldExit())
        {
            sort (list);
            completionCallback (list);
        }

        return jobHasFinished;
    }

    Array<File> pathsToScan;
    std::function<void (const ModuleIDAndFolderList&)> completionCallback;
};

AvailableModuleList::AvailableModuleList()
{
}

ThreadPoolJob* AvailableModuleList::createScannerJob (const Array<File>& paths)
{
    return new ModuleScannerJob (paths, [this] (ModuleIDAndFolderList scannedModuleList)
                                        {
                                            {
                                                const ScopedLock swapLock (lock);
                                                moduleList.swap (scannedModuleList);
                                            }

                                            listeners.call ([] (Listener& l) { MessageManager::callAsync ([&] { l.availableModulesChanged(); }); });
                                        });
}

void AvailableModuleList::removePendingAndAddJob (ThreadPoolJob* jobToAdd)
{
    scanPool.removeAllJobs (false, 100);
    scanPool.addJob (jobToAdd, true);
}

void AvailableModuleList::scanPaths (const Array<File>& paths)
{
    auto* job = createScannerJob (paths);

    removePendingAndAddJob (job);
    scanPool.waitForJobToFinish (job, -1);
}

void AvailableModuleList::scanPathsAsync (const Array<File>& paths)
{
    removePendingAndAddJob (createScannerJob (paths));
}

ModuleIDAndFolderList AvailableModuleList::getAllModules() const
{
    const ScopedLock readLock (lock);

    return moduleList;
}

ModuleIDAndFolder AvailableModuleList::getModuleWithID (const String& id) const
{
    const ScopedLock readLock (lock);

    for (auto& mod : moduleList)
        if (mod.first == id)
            return mod;

    return {};
}

void AvailableModuleList::removeDuplicates (const ModuleIDAndFolderList& other)
{
    const ScopedLock readLock (lock);

    for (auto& m : other)
    {
        auto pos = std::find (moduleList.begin(), moduleList.end(), m);

        if (pos != moduleList.end())
            moduleList.erase (pos);
    }
}

//==============================================================================
LibraryModule::LibraryModule (const ModuleDescription& d)
    : moduleInfo (d)
{
}

//==============================================================================
void LibraryModule::writeIncludes (ProjectSaver& projectSaver, OutputStream& out)
{
    auto& project = projectSaver.project;
    auto& modules = project.getEnabledModules();

    auto id = getID();

    if (modules.shouldCopyModuleFilesLocally (id).getValue())
    {
        auto juceModuleFolder = moduleInfo.getFolder();

        auto localModuleFolder = project.getLocalModuleFolder (id);
        localModuleFolder.createDirectory();
        projectSaver.copyFolder (juceModuleFolder, localModuleFolder);
    }

    out << "#include <" << moduleInfo.moduleFolder.getFileName() << "/"
        << moduleInfo.getHeader().getFileName()
        << ">" << newLine;
}

//==============================================================================
static void parseAndAddLibs (StringArray& libList, const String& libs)
{
    libList.addTokens (libs, ", ", {});
    libList.trim();
    libList.sort (false);
    libList.removeDuplicates (false);
}

void LibraryModule::addSettingsForModuleToExporter (ProjectExporter& exporter, ProjectSaver& projectSaver) const
{
    auto& project = exporter.getProject();

    auto moduleRelativePath = exporter.getModuleFolderRelativeToProject (getID());

    exporter.addToExtraSearchPaths (moduleRelativePath.getParentDirectory());

    String libDirPlatform;
    if (exporter.isLinux())
        libDirPlatform = "Linux";
    else if (exporter.isCodeBlocks() && exporter.isWindows())
        libDirPlatform = "MinGW";
    else
        libDirPlatform = exporter.getTargetFolder().getFileName();

    auto libSubdirPath = moduleRelativePath.toUnixStyle() + "/libs/" + libDirPlatform;
    auto moduleLibDir = File (project.getProjectFolder().getFullPathName() + "/" + libSubdirPath);

    if (moduleLibDir.exists())
        exporter.addToModuleLibPaths ({ libSubdirPath, moduleRelativePath.getRoot() });

    auto extraInternalSearchPaths = moduleInfo.getExtraSearchPaths().trim();

    if (extraInternalSearchPaths.isNotEmpty())
    {
        auto paths = StringArray::fromTokens (extraInternalSearchPaths, true);

        for (auto& path : paths)
            exporter.addToExtraSearchPaths (moduleRelativePath.getChildFile (path.unquoted()));
    }

    {
        auto extraDefs = moduleInfo.getPreprocessorDefs().trim();

        if (extraDefs.isNotEmpty())
            exporter.getExporterPreprocessorDefsValue() = exporter.getExporterPreprocessorDefsString() + "\n" + extraDefs;
    }

    {
        Array<File> compiled;
        auto& modules = project.getEnabledModules();
        auto id = getID();

        auto localModuleFolder = modules.shouldCopyModuleFilesLocally (id).getValue() ? project.getLocalModuleFolder (id)
                                                                                      : moduleInfo.getFolder();

        findAndAddCompiledUnits (exporter, &projectSaver, compiled);

        if (modules.shouldShowAllModuleFilesInProject (id).getValue())
            addBrowseableCode (exporter, compiled, localModuleFolder);
    }

    if (exporter.isXcode())
    {
        auto& xcodeExporter = dynamic_cast<XcodeProjectExporter&> (exporter);

        if (project.isAUPluginHost())
            xcodeExporter.xcodeFrameworks.addTokens (xcodeExporter.isOSX() ? "AudioUnit CoreAudioKit" : "CoreAudioKit", false);

        auto frameworks = moduleInfo.moduleInfo [xcodeExporter.isOSX() ? "OSXFrameworks" : "iOSFrameworks"].toString();
        xcodeExporter.xcodeFrameworks.addTokens (frameworks, ", ", {});

        parseAndAddLibs (xcodeExporter.xcodeLibs, moduleInfo.moduleInfo [exporter.isOSX() ? "OSXLibs" : "iOSLibs"].toString());
    }
    else if (exporter.isLinux())
    {
        parseAndAddLibs (exporter.linuxLibs, moduleInfo.moduleInfo ["linuxLibs"].toString());
        parseAndAddLibs (exporter.linuxPackages, moduleInfo.moduleInfo ["linuxPackages"].toString());
    }
    else if (exporter.isWindows())
    {
        if (exporter.isCodeBlocks())
            parseAndAddLibs (exporter.mingwLibs, moduleInfo.moduleInfo ["mingwLibs"].toString());
        else
            parseAndAddLibs (exporter.windowsLibs, moduleInfo.moduleInfo ["windowsLibs"].toString());
    }
    else if (exporter.isAndroid())
    {
        parseAndAddLibs (exporter.androidLibs, moduleInfo.moduleInfo ["androidLibs"].toString());
    }
}

void LibraryModule::getConfigFlags (Project& project, OwnedArray<Project::ConfigFlag>& flags) const
{
    auto header = moduleInfo.getHeader();
    jassert (header.exists());

    StringArray lines;
    header.readLines (lines);

    for (int i = 0; i < lines.size(); ++i)
    {
        auto line = lines[i].trim();

        if (line.startsWith ("/**") && line.containsIgnoreCase ("Config:"))
        {
            std::unique_ptr<Project::ConfigFlag> config (new Project::ConfigFlag());
            config->sourceModuleID = getID();
            config->symbol = line.fromFirstOccurrenceOf (":", false, false).trim();

            if (config->symbol.length() > 2)
            {
                ++i;

                while (! (lines[i].contains ("*/") || lines[i].contains ("@see")))
                {
                    if (lines[i].trim().isNotEmpty())
                        config->description = config->description.trim() + " " + lines[i].trim();

                    ++i;
                }

                config->description = config->description.upToFirstOccurrenceOf ("*/", false, false);
                config->value = project.getConfigFlag (config->symbol);

                i += 2;
                if (lines[i].contains ("#define " + config->symbol))
                {
                    auto value = lines[i].fromFirstOccurrenceOf ("#define " + config->symbol, false, true).trim();
                    config->value.setDefault (value == "0" ? false : true);
                }

                auto currentValue = config->value.get().toString();

                if      (currentValue == "enabled")     config->value = true;
                else if (currentValue == "disabled")    config->value = false;

                flags.add (config.release());
            }
        }
    }
}

//==============================================================================
struct FileSorter
{
    static int compareElements (const File& f1, const File& f2)
    {
        return f1.getFileName().compareNatural (f2.getFileName());
    }
};

bool LibraryModule::CompileUnit::hasSuffix (const File& f, const char* suffix)
{
    auto fileWithoutSuffix = f.getFileNameWithoutExtension() + ".";

    return fileWithoutSuffix.containsIgnoreCase (suffix + String ("."))
             || fileWithoutSuffix.containsIgnoreCase (suffix + String ("_"));
}

void LibraryModule::CompileUnit::writeInclude (MemoryOutputStream&) const
{
}

bool LibraryModule::CompileUnit::isNeededForExporter (ProjectExporter& exporter) const
{
    if ((hasSuffix (file, "_OSX")        && ! exporter.isOSX())
     || (hasSuffix (file, "_iOS")        && ! exporter.isiOS())
     || (hasSuffix (file, "_Windows")    && ! exporter.isWindows())
     || (hasSuffix (file, "_Linux")      && ! exporter.isLinux())
     || (hasSuffix (file, "_Android")    && ! exporter.isAndroid()))
        return false;

    auto targetType = Project::getTargetTypeFromFilePath (file, false);

    if (targetType != ProjectType::Target::unspecified && ! exporter.shouldBuildTargetType (targetType))
        return false;

    return exporter.usesMMFiles() ? isCompiledForObjC
                                  : isCompiledForNonObjC;
}

String LibraryModule::CompileUnit::getFilenameForProxyFile() const
{
    return "include_" + file.getFileName();
}

Array<LibraryModule::CompileUnit> LibraryModule::getAllCompileUnits (ProjectType::Target::Type forTarget) const
{
    auto files = getFolder().findChildFiles (File::findFiles, false);

    FileSorter sorter;
    files.sort (sorter);

    Array<LibraryModule::CompileUnit> units;

    for (auto& file : files)
    {
        if (file.getFileName().startsWithIgnoreCase (getID())
              && file.hasFileExtension (sourceFileExtensions))
        {
            if (forTarget == ProjectType::Target::unspecified
             || forTarget == Project::getTargetTypeFromFilePath (file, true))
            {
                CompileUnit cu;
                cu.file = file;
                units.add (cu);
            }
        }
    }

    for (auto& cu : units)
    {
        cu.isCompiledForObjC = true;
        cu.isCompiledForNonObjC = ! cu.file.hasFileExtension ("mm;m");

        if (cu.isCompiledForNonObjC)
            if (files.contains (cu.file.withFileExtension ("mm")))
                cu.isCompiledForObjC = false;

        jassert (cu.isCompiledForObjC || cu.isCompiledForNonObjC);
    }

    return units;
}

void LibraryModule::findAndAddCompiledUnits (ProjectExporter& exporter,
                                             ProjectSaver* projectSaver,
                                             Array<File>& result,
                                             ProjectType::Target::Type forTarget) const
{
    for (auto& cu : getAllCompileUnits (forTarget))
    {
        if (cu.isNeededForExporter (exporter))
        {
            auto localFile = exporter.getProject().getGeneratedCodeFolder()
                                                  .getChildFile (cu.getFilenameForProxyFile());
            result.add (localFile);

            if (projectSaver != nullptr)
                projectSaver->addFileToGeneratedGroup (localFile);
        }
    }
}

static void addFileWithGroups (Project::Item& group, const RelativePath& file, const String& path)
{
    auto slash = path.indexOfChar (File::getSeparatorChar());

    if (slash >= 0)
    {
        auto topLevelGroup = path.substring (0, slash);
        auto remainingPath = path.substring (slash + 1);

        auto newGroup = group.getOrCreateSubGroup (topLevelGroup);
        addFileWithGroups (newGroup, file, remainingPath);
    }
    else
    {
        if (! group.containsChildForFile (file))
            group.addRelativeFile (file, -1, false);
    }
}

void LibraryModule::findBrowseableFiles (const File& folder, Array<File>& filesFound) const
{
    Array<File> tempList;
    FileSorter sorter;

    DirectoryIterator iter (folder, true, "*", File::findFiles);
    bool isHiddenFile;

    while (iter.next (nullptr, &isHiddenFile, nullptr, nullptr, nullptr, nullptr))
        if (! isHiddenFile && iter.getFile().hasFileExtension (browseableFileExtensions))
            tempList.addSorted (sorter, iter.getFile());

    filesFound.addArray (tempList);
}

void LibraryModule::addBrowseableCode (ProjectExporter& exporter, const Array<File>& compiled, const File& localModuleFolder) const
{
    if (sourceFiles.isEmpty())
        findBrowseableFiles (localModuleFolder, sourceFiles);

    auto sourceGroup = Project::Item::createGroup (exporter.getProject(), getID(), "__mainsourcegroup" + getID(), false);

    auto moduleFromProject = exporter.getModuleFolderRelativeToProject (getID());

    auto moduleHeader = moduleInfo.getHeader();

    for (auto& sourceFile : sourceFiles)
    {
        auto pathWithinModule = FileHelpers::getRelativePathFrom (sourceFile, localModuleFolder);

        // (Note: in exporters like MSVC we have to avoid adding the same file twice, even if one of those instances
        // is flagged as being excluded from the build, because this overrides the other and it fails to compile)
        if ((exporter.canCopeWithDuplicateFiles() || ! compiled.contains (sourceFile)) && sourceFile != moduleHeader)
            addFileWithGroups (sourceGroup,
                               moduleFromProject.getChildFile (pathWithinModule),
                               pathWithinModule);
    }

    sourceGroup.sortAlphabetically (true, true);
    sourceGroup.addFileAtIndex (moduleHeader, -1, false);

    exporter.getModulesGroup().state.appendChild (sourceGroup.state.createCopy(), nullptr);
}


//==============================================================================
EnabledModuleList::EnabledModuleList (Project& p, const ValueTree& s)
    : project (p), state (s)
{
}

ModuleDescription EnabledModuleList::getModuleInfo (const String& moduleID)
{
    return ModuleDescription (project.getModuleWithID (moduleID).second);
}

bool EnabledModuleList::isModuleEnabled (const String& moduleID) const
{
    return state.getChildWithProperty (Ids::ID, moduleID).isValid();
}

bool EnabledModuleList::isAudioPluginModuleMissing() const
{
    return project.getProjectType().isAudioPlugin()
            && ! isModuleEnabled ("juce_audio_plugin_client");
}

bool EnabledModuleList::shouldUseGlobalPath (const String& moduleID) const
{
    return static_cast<bool> (state.getChildWithProperty (Ids::ID, moduleID)
                                   .getProperty (Ids::useGlobalPath));
}

Value EnabledModuleList::getShouldUseGlobalPathValue (const String& moduleID) const
{
    return state.getChildWithProperty (Ids::ID, moduleID)
                .getPropertyAsValue (Ids::useGlobalPath, getUndoManager());
}

Value EnabledModuleList::shouldShowAllModuleFilesInProject (const String& moduleID)
{
    return state.getChildWithProperty (Ids::ID, moduleID)
                .getPropertyAsValue (Ids::showAllCode, getUndoManager());
}

struct ModuleTreeSorter
{
    static int compareElements (const ValueTree& m1, const ValueTree& m2)
    {
        return m1[Ids::ID].toString().compareIgnoreCase (m2[Ids::ID]);
    }
};

void EnabledModuleList::sortAlphabetically()
{
    ModuleTreeSorter sorter;
    state.sort (sorter, getUndoManager(), false);
}

Value EnabledModuleList::shouldCopyModuleFilesLocally (const String& moduleID) const
{
    return state.getChildWithProperty (Ids::ID, moduleID)
                .getPropertyAsValue (Ids::useLocalCopy, getUndoManager());
}

void EnabledModuleList::addModule (const File& moduleFolder, bool copyLocally, bool useGlobalPath, bool sendAnalyticsEvent)
{
    ModuleDescription info (moduleFolder);

    if (info.isValid())
    {
        auto moduleID = info.getID();

        if (! isModuleEnabled (moduleID))
        {
            ValueTree module (Ids::MODULE);
            module.setProperty (Ids::ID, moduleID, getUndoManager());

            state.appendChild (module, getUndoManager());
            sortAlphabetically();

            shouldShowAllModuleFilesInProject (moduleID) = true;
            shouldCopyModuleFilesLocally (moduleID) = copyLocally;
            getShouldUseGlobalPathValue (moduleID) = useGlobalPath;

            RelativePath path (moduleFolder.getParentDirectory(),
                               project.getProjectFolder(), RelativePath::projectFolder);

            for (Project::ExporterIterator exporter (project); exporter.next();)
                exporter->getPathForModuleValue (moduleID) = path.toUnixStyle();

            if (! useGlobalPath)
                project.rescanExporterPathModules (false);

            if (sendAnalyticsEvent)
            {
                StringPairArray data;
                data.set ("label", moduleID);

                Analytics::getInstance()->logEvent ("Module Added", data, ProjucerAnalyticsEvent::projectEvent);
            }
        }
    }
}

void EnabledModuleList::removeModule (String moduleID) // must be pass-by-value, and not a const ref!
{
    for (auto i = state.getNumChildren(); --i >= 0;)
        if (state.getChild(i) [Ids::ID] == moduleID)
            state.removeChild (i, getUndoManager());

    for (Project::ExporterIterator exporter (project); exporter.next();)
        exporter->removePathForModule (moduleID);
}

void EnabledModuleList::createRequiredModules (OwnedArray<LibraryModule>& modules)
{
    for (int i = 0; i < getNumModules(); ++i)
        modules.add (new LibraryModule (getModuleInfo (getModuleID (i))));
}

StringArray EnabledModuleList::getAllModules() const
{
    StringArray moduleIDs;

    for (int i = 0; i < getNumModules(); ++i)
        moduleIDs.add (getModuleID (i));

    return moduleIDs;
}

static void getDependencies (Project& project, const String& moduleID, StringArray& dependencies)
{
    auto info = project.getEnabledModules().getModuleInfo (moduleID);

    for (auto uid : info.getDependencies())
    {
        if (! dependencies.contains (uid, true))
        {
            dependencies.add (uid);
            getDependencies (project, uid, dependencies);
        }
    }
}

StringArray EnabledModuleList::getExtraDependenciesNeeded (const String& moduleID) const
{
    StringArray dependencies, extraDepsNeeded;
    getDependencies (project, moduleID, dependencies);

    for (auto dep : dependencies)
        if (dep != moduleID && ! isModuleEnabled (dep))
            extraDepsNeeded.add (dep);

    return extraDepsNeeded;
}

bool EnabledModuleList::doesModuleHaveHigherCppStandardThanProject (const String& moduleID)
{
    auto projectCppStandard = project.getCppStandardString();

    if (projectCppStandard == "latest")
        return false;

    auto moduleCppStandard = getModuleInfo (moduleID).getMinimumCppStandard();

    return (moduleCppStandard.getIntValue() > projectCppStandard.getIntValue());
}

bool EnabledModuleList::areMostModulesUsingGlobalPath() const
{
    int numYes = 0, numNo = 0;

    for (auto i = getNumModules(); --i >= 0;)
    {
        if (shouldUseGlobalPath (getModuleID (i)))
            ++numYes;
        else
            ++numNo;
    }

    return numYes > numNo;
}

bool EnabledModuleList::areMostModulesCopiedLocally() const
{
    int numYes = 0, numNo = 0;

    for (auto i = getNumModules(); --i >= 0;)
    {
        if (shouldCopyModuleFilesLocally (getModuleID (i)).getValue())
            ++numYes;
        else
            ++numNo;
    }

    return numYes > numNo;
}

void EnabledModuleList::setLocalCopyModeForAllModules (bool copyLocally)
{
    for (auto i = getNumModules(); --i >= 0;)
        shouldCopyModuleFilesLocally (project.getEnabledModules().getModuleID (i)) = copyLocally;
}

File EnabledModuleList::findDefaultModulesFolder (Project& project)
{
    File globalPath (getAppSettings().getStoredPath (Ids::defaultJuceModulePath, TargetOS::getThisOS()).get().toString());

    if (globalPath.exists())
        return globalPath;

    for (auto& exporterPathModule : project.getExporterPathsModuleList().getAllModules())
    {
        auto f = exporterPathModule.second;

        if (f.isDirectory())
            return f.getParentDirectory();
    }

    return File::getCurrentWorkingDirectory();
}

void EnabledModuleList::addModuleFromUserSelectedFile()
{
    static auto lastLocation = findDefaultModulesFolder (project);

    FileChooser fc ("Select a module to add...", lastLocation, {});

    if (fc.browseForDirectory())
    {
        lastLocation = fc.getResult();
        addModuleOfferingToCopy (lastLocation, true);
    }
}

void EnabledModuleList::addModuleInteractive (const String& moduleID)
{
    auto f = project.getModuleWithID (moduleID).second;

    if (f != File())
    {
        addModule (f, areMostModulesCopiedLocally(), areMostModulesUsingGlobalPath(), true);
        return;
    }

    addModuleFromUserSelectedFile();
}

void EnabledModuleList::addModuleOfferingToCopy (const File& f, bool isFromUserSpecifiedFolder)
{
    ModuleDescription m (f);

    if (! m.isValid())
    {
        AlertWindow::showMessageBoxAsync (AlertWindow::InfoIcon,
                                          "Add Module", "This wasn't a valid module folder!");
        return;
    }

    if (isModuleEnabled (m.getID()))
    {
        AlertWindow::showMessageBoxAsync (AlertWindow::InfoIcon,
                                          "Add Module", "The project already contains this module!");
        return;
    }

    addModule (m.moduleFolder, areMostModulesCopiedLocally(),
               isFromUserSpecifiedFolder ? false : areMostModulesUsingGlobalPath(),
               true);
}

bool isJUCEFolder (const File& f)
{
    return isJUCEModulesFolder (f.getChildFile ("modules"));
}

bool isJUCEModulesFolder (const File& f)
{
    return f.isDirectory() && f.getChildFile ("juce_core").isDirectory();
}
