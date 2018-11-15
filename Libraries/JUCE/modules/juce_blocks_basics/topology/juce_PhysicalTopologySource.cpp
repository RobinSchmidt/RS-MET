/*
  ==============================================================================

   This file is part of the JUCE library.
   Copyright (c) 2017 - ROLI Ltd.

   JUCE is an open source library subject to commercial or open-source
   licensing.

   The code included in this file is provided under the terms of the ISC license
   http://www.isc.org/downloads/software-support-policy/isc-license. Permission
   To use, copy, modify, and/or distribute this software for any purpose with or
   without fee is hereby granted provided that the above copyright notice and
   this permission notice appear in all copies.

   JUCE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY, AND ALL WARRANTIES, WHETHER
   EXPRESSED OR IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR PURPOSE, ARE
   DISCLAIMED.

  ==============================================================================
*/

namespace juce
{

//==============================================================================
/** These can be useful when debugging the topology. */
#define LOG_BLOCKS_CONNECTIVITY 0
#define LOG_BLOCKS_PINGS 0
#define DUMP_BANDWIDTH_STATS 0

#define TOPOLOGY_LOG(text) \
    JUCE_BLOCK_WITH_FORCED_SEMICOLON (juce::String buf ("Topology Src:   "); \
    juce::Logger::outputDebugString (buf << text);)

#if LOG_BLOCKS_CONNECTIVITY
 #define LOG_CONNECTIVITY(text) TOPOLOGY_LOG(text)
#else
 #define LOG_CONNECTIVITY(text)
#endif

#if LOG_BLOCKS_PINGS
 #define LOG_PING(text) TOPOLOGY_LOG(text)
#else
 #define LOG_PING(text)
#endif

#if DUMP_BANDWIDTH_STATS
namespace
{
    struct PortIOStats
    {
        PortIOStats (const char* nm) : name (nm) {}

        const char* const name;
        int byteCount           = 0;
        int messageCount        = 0;
        int bytesPerSec         = 0;
        int largestMessageBytes = 0;
        int lastMessageBytes    = 0;

        void update (double elapsedSec)
        {
            if (byteCount > 0)
            {
                bytesPerSec = (int) (byteCount / elapsedSec);
                byteCount = 0;
                juce::Logger::writeToLog (getString());
            }
        }

        juce::String getString() const
        {
            return juce::String (name) + ": "
                    + "count="    + juce::String (messageCount).paddedRight (' ', 7)
                    + "rate="     + (juce::String (bytesPerSec / 1024.0f, 1) + " Kb/sec").paddedRight (' ', 11)
                    + "largest="  + (juce::String (largestMessageBytes) + " bytes").paddedRight (' ', 11)
                    + "last="     + (juce::String (lastMessageBytes) + " bytes").paddedRight (' ', 11);
        }

        void registerMessage (int numBytes) noexcept
        {
            byteCount += numBytes;
            ++messageCount;
            lastMessageBytes = numBytes;
            largestMessageBytes = juce::jmax (largestMessageBytes, numBytes);
        }
    };

    static PortIOStats inputStats { "Input" }, outputStats { "Output" };
    static uint32 startTime = 0;

    static inline void resetOnSecondBoundary()
    {
        auto now = juce::Time::getMillisecondCounter();
        double elapsedSec = (now - startTime) / 1000.0;

        if (elapsedSec >= 1.0)
        {
            inputStats.update (elapsedSec);
            outputStats.update (elapsedSec);
            startTime = now;
        }
    }

    static inline void registerBytesOut (int numBytes)
    {
        outputStats.registerMessage (numBytes);
        resetOnSecondBoundary();
    }

    static inline void registerBytesIn (int numBytes)
    {
        inputStats.registerMessage (numBytes);
        resetOnSecondBoundary();
    }
}

juce::String getMidiIOStats()
{
    return inputStats.getString() + "   " + outputStats.getString();
}
#endif

//==============================================================================
struct PhysicalTopologySource::Internal
{
    struct Detector;
    struct BlockImplementation;
    struct ControlButtonImplementation;
    struct RotaryDialImplementation;
    struct TouchSurfaceImplementation;
    struct LEDGridImplementation;
    struct LEDRowImplementation;

    //==============================================================================
    struct MIDIDeviceConnection  : public DeviceConnection,
                                   public juce::MidiInputCallback
    {
        MIDIDeviceConnection() {}

        ~MIDIDeviceConnection()
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

            listeners.call ([this] (Listener& l) { l.connectionBeingDeleted (*this); });

            if (midiInput != nullptr)
                midiInput->stop();

            if (interprocessLock != nullptr)
                interprocessLock->exit();
        }

        bool lockAgainstOtherProcesses (const String& midiInName, const String& midiOutName)
        {
            interprocessLock.reset (new juce::InterProcessLock ("blocks_sdk_"
                                                                  + File::createLegalFileName (midiInName)
                                                                  + "_" + File::createLegalFileName (midiOutName)));
            if (interprocessLock->enter (500))
                return true;

            interprocessLock = nullptr;
            return false;
        }

        struct Listener
        {
            virtual ~Listener() {}

            virtual void handleIncomingMidiMessage (const juce::MidiMessage& message) = 0;
            virtual void connectionBeingDeleted (const MIDIDeviceConnection&) = 0;
        };

        void addListener (Listener* l)
        {
            listeners.add (l);
        }

        void removeListener (Listener* l)
        {
            listeners.remove (l);
        }

        bool sendMessageToDevice (const void* data, size_t dataSize) override
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED // This method must only be called from the message thread!

            jassert (dataSize > sizeof (BlocksProtocol::roliSysexHeader) + 2);
            jassert (memcmp (data, BlocksProtocol::roliSysexHeader, sizeof (BlocksProtocol::roliSysexHeader)) == 0);
            jassert (static_cast<const uint8*> (data)[dataSize - 1] == 0xf7);

            if (midiOutput != nullptr)
            {
                midiOutput->sendMessageNow (juce::MidiMessage (data, (int) dataSize));
                return true;
            }

            return false;
        }

        void handleIncomingMidiMessage (juce::MidiInput*, const juce::MidiMessage& message) override
        {
            const auto data = message.getRawData();
            const int dataSize = message.getRawDataSize();
            const int bodySize = dataSize - (int) (sizeof (BlocksProtocol::roliSysexHeader) + 1);

            if (bodySize > 0 && memcmp (data, BlocksProtocol::roliSysexHeader, sizeof (BlocksProtocol::roliSysexHeader)) == 0)
                if (handleMessageFromDevice != nullptr)
                    handleMessageFromDevice (data + sizeof (BlocksProtocol::roliSysexHeader), (size_t) bodySize);

            listeners.call ([&] (Listener& l) { l.handleIncomingMidiMessage (message); });
        }

        std::unique_ptr<juce::MidiInput> midiInput;
        std::unique_ptr<juce::MidiOutput> midiOutput;

    private:
        juce::ListenerList<Listener> listeners;
        std::unique_ptr<juce::InterProcessLock> interprocessLock;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MIDIDeviceConnection)
    };

    //==============================================================================
    struct MIDIDeviceDetector  : public DeviceDetector
    {
        MIDIDeviceDetector() {}

        juce::StringArray scanForDevices() override
        {
            juce::StringArray result;

            for (auto& pair : findDevices())
                result.add (pair.inputName + " & " + pair.outputName);

            return result;
        }

        DeviceConnection* openDevice (int index) override
        {
            auto pair = findDevices()[index];

            if (pair.inputIndex >= 0 && pair.outputIndex >= 0)
            {
                std::unique_ptr<MIDIDeviceConnection> dev (new MIDIDeviceConnection());

                if (dev->lockAgainstOtherProcesses (pair.inputName, pair.outputName))
                {
                    lockedFromOutside = false;

                    dev->midiInput.reset (juce::MidiInput::openDevice (pair.inputIndex, dev.get()));
                    dev->midiOutput.reset (juce::MidiOutput::openDevice (pair.outputIndex));

                    if (dev->midiInput != nullptr)
                    {
                        dev->midiInput->start();
                        return dev.release();
                    }
                }
                else
                {
                    lockedFromOutside = true;
                }
            }

            return nullptr;
        }

        bool isLockedFromOutside() const override
        {
            return lockedFromOutside && ! findDevices().isEmpty();
        }

        static bool isBlocksMidiDeviceName (const juce::String& name)
        {
            return name.indexOf (" BLOCK") > 0 || name.indexOf (" Block") > 0;
        }

        static String cleanBlocksDeviceName (juce::String name)
        {
            name = name.trim();

            if (name.endsWith (" IN)"))
                return name.dropLastCharacters (4);

            if (name.endsWith (" OUT)"))
                return name.dropLastCharacters (5);

            const int openBracketPosition = name.lastIndexOfChar ('[');
            if (openBracketPosition != -1 && name.endsWith ("]"))
                return name.dropLastCharacters (name.length() - openBracketPosition);

            return name;
        }

        struct MidiInputOutputPair
        {
            juce::String outputName, inputName;
            int outputIndex = -1, inputIndex = -1;
        };

        static juce::Array<MidiInputOutputPair> findDevices()
        {
            juce::Array<MidiInputOutputPair> result;

            auto midiInputs  = juce::MidiInput::getDevices();
            auto midiOutputs = juce::MidiOutput::getDevices();

            for (int j = 0; j < midiInputs.size(); ++j)
            {
                if (isBlocksMidiDeviceName (midiInputs[j]))
                {
                    MidiInputOutputPair pair;
                    pair.inputName = midiInputs[j];
                    pair.inputIndex = j;

                    String cleanedInputName = cleanBlocksDeviceName (pair.inputName);
                    for (int i = 0; i < midiOutputs.size(); ++i)
                    {
                        if (cleanBlocksDeviceName (midiOutputs[i]) == cleanedInputName)
                        {
                            pair.outputName = midiOutputs[i];
                            pair.outputIndex = i;
                            break;
                        }
                    }

                    result.add (pair);
                }
            }

            return result;
        }

    private:
        bool lockedFromOutside = true;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MIDIDeviceDetector)
    };

    //==============================================================================
    struct DeviceInfo
    {
        // VS2015 requires a constructor to avoid aggregate initialization
        DeviceInfo (Block::UID buid, BlocksProtocol::TopologyIndex tidx, BlocksProtocol::BlockSerialNumber s,
                    BlocksProtocol::VersionNumber v, BlocksProtocol::BlockName n, bool master = false)
            : uid (buid), index (tidx), serial (s), version (v), name (n), isMaster (master)
        {
        }

        Block::UID uid {};
        BlocksProtocol::TopologyIndex index;
        BlocksProtocol::BlockSerialNumber serial;
        BlocksProtocol::VersionNumber version;
        BlocksProtocol::BlockName name;
        bool isMaster {};
    };

    static juce::String getVersionString (const BlocksProtocol::VersionNumber& v)
    {
        return juce::String (reinterpret_cast<const char*> (v.version),
                             std::min (sizeof (v.version), static_cast<size_t> (v.length)));
    }

    static juce::String getNameString (const BlocksProtocol::BlockName& n)
    {
        return juce::String (reinterpret_cast<const char*> (n.name),
                             std::min (sizeof (n.name), static_cast<size_t> (n.length)));
    }

    static Block::Timestamp deviceTimestampToHost (uint32 timestamp) noexcept
    {
        return static_cast<Block::Timestamp> (timestamp);
    }

    static juce::Array<DeviceInfo> getArrayOfDeviceInfo (const juce::Array<BlocksProtocol::DeviceStatus>& devices)
    {
        juce::Array<DeviceInfo> result;
        bool isFirst = true; // TODO: First block not always master block! Assumption violated.

        for (auto& device : devices)
        {
            BlocksProtocol::VersionNumber version;
            BlocksProtocol::BlockName name;

            result.add ({ getBlockUIDFromSerialNumber (device.serialNumber),
                          device.index,
                          device.serialNumber,
                          version,
                          name,
                          isFirst });

            isFirst = false;
        }

        return result;
    }

    static bool containsBlockWithUID (const Block::Array& blocks, Block::UID uid) noexcept
    {
        for (auto&& block : blocks)
            if (block->uid == uid)
                return true;

        return false;
    }

    static bool versionNumberChanged (const DeviceInfo& device, juce::String version) noexcept
    {
        auto deviceVersion = getVersionString (device.version);
        return deviceVersion != version && deviceVersion.isNotEmpty();
    }

    static bool nameIsValid (const DeviceInfo& device)
    {
        return device.name.length > 0;
    }

    static void setVersionNumberForBlock (const DeviceInfo& deviceInfo, Block& block) noexcept
    {
        jassert (deviceInfo.uid == block.uid);
        block.versionNumber = getVersionString (deviceInfo.version);
    }

    static void setNameForBlock (const DeviceInfo& deviceInfo, Block& block)
    {
        jassert (deviceInfo.uid == block.uid);
        block.name = getNameString (deviceInfo.name);
    }

    //==============================================================================
    struct ConnectedDeviceGroup  : private juce::AsyncUpdater,
                                   private juce::Timer
    {
        ConnectedDeviceGroup (Detector& d, const juce::String& name, DeviceConnection* connection)
            : detector (d), deviceName (name), deviceConnection (connection)
        {
            deviceConnection->handleMessageFromDevice = [this] (const void* data, size_t dataSize)
            {
                this->handleIncomingMessage (data, dataSize);
            };

            startTimer (200);
            sendTopologyRequest();
        }

        bool isStillConnected (const juce::StringArray& detectedDevices) const noexcept
        {
            return detectedDevices.contains (deviceName)
                && ! failedToGetTopology();
        }

        int getIndexFromDeviceID (Block::UID uid) const noexcept
        {
            for (auto& d : currentDeviceInfo)
                if (d.uid == uid)
                    return d.index;

            return -1;
        }

        const DeviceInfo* getDeviceInfoFromUID (Block::UID uid) const noexcept
        {
            for (auto& d : currentDeviceInfo)
                if (d.uid == uid)
                    return &d;

            return nullptr;
        }

        const BlocksProtocol::DeviceStatus* getLastStatus (Block::UID deviceID) const noexcept
        {
            for (auto&& status : currentTopologyDevices)
                if (getBlockUIDFromSerialNumber (status.serialNumber) == deviceID)
                    return &status;

            return nullptr;
        }

        void notifyBlockIsRestarting (Block::UID deviceID)
        {
            forceApiDisconnected (deviceID);
        }

        //==============================================================================
        // The following methods will be called by the HostPacketDecoder:
        void beginTopology (int numDevices, int numConnections)
        {
            incomingTopologyDevices.clearQuick();
            incomingTopologyDevices.ensureStorageAllocated (numDevices);
            incomingTopologyConnections.clearQuick();
            incomingTopologyConnections.ensureStorageAllocated (numConnections);
        }

        void extendTopology (int numDevices, int numConnections)
        {
            incomingTopologyDevices.ensureStorageAllocated (incomingTopologyDevices.size() + numDevices);
            incomingTopologyConnections.ensureStorageAllocated (incomingTopologyConnections.size() + numConnections);
        }

        void handleTopologyDevice (BlocksProtocol::DeviceStatus status)
        {
            incomingTopologyDevices.add (status);
        }

        void handleTopologyConnection (BlocksProtocol::DeviceConnection connection)
        {
            incomingTopologyConnections.add (connection);
        }

        void endTopology()
        {
            currentDeviceInfo = getArrayOfDeviceInfo (incomingTopologyDevices);
            currentDeviceConnections = getArrayOfConnections (incomingTopologyConnections);
            currentTopologyDevices = incomingTopologyDevices;
            lastTopologyReceiveTime = juce::Time::getCurrentTime();

            const int numRemoved = blockPings.removeIf ([this] (auto& ping)
            {
                for (auto& info : currentDeviceInfo)
                    if (info.uid == ping.blockUID)
                        return false;

                LOG_CONNECTIVITY ("API Disconnected by topology update " << ping.blockUID);
                return true;
            });

            if (numRemoved > 0)
                detector.handleTopologyChange();
        }

        void handleVersion (BlocksProtocol::DeviceVersion version)
        {
            for (auto& d : currentDeviceInfo)
                if (d.index == version.index && version.version.length > 1)
                    d.version = version.version;
        }

        void handleName (BlocksProtocol::DeviceName name)
        {
            for (auto& d : currentDeviceInfo)
                if (d.index == name.index && name.name.length > 1)
                    d.name = name.name;
        }

        void handleControlButtonUpDown (BlocksProtocol::TopologyIndex deviceIndex, uint32 timestamp,
                                        BlocksProtocol::ControlButtonID buttonID, bool isDown)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleButtonChange (deviceID, deviceTimestampToHost (timestamp), buttonID.get(), isDown);
        }

        void handleCustomMessage (BlocksProtocol::TopologyIndex deviceIndex, uint32 timestamp, const int32* data)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleCustomMessage (deviceID, deviceTimestampToHost (timestamp), data);
        }

        void handleTouchChange (BlocksProtocol::TopologyIndex deviceIndex,
                                uint32 timestamp,
                                BlocksProtocol::TouchIndex touchIndex,
                                BlocksProtocol::TouchPosition position,
                                BlocksProtocol::TouchVelocity velocity,
                                bool isStart, bool isEnd)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
            {
                TouchSurface::Touch touch;

                touch.index             = (int) touchIndex.get();
                touch.x                 = position.x.toUnipolarFloat();
                touch.y                 = position.y.toUnipolarFloat();
                touch.z                 = position.z.toUnipolarFloat();
                touch.xVelocity         = velocity.vx.toBipolarFloat();
                touch.yVelocity         = velocity.vy.toBipolarFloat();
                touch.zVelocity         = velocity.vz.toBipolarFloat();
                touch.eventTimestamp    = deviceTimestampToHost (timestamp);
                touch.isTouchStart      = isStart;
                touch.isTouchEnd        = isEnd;
                touch.blockUID          = deviceID;

                setTouchStartPosition (touch);

                detector.handleTouchChange (deviceID, touch);
            }
        }

        void setTouchStartPosition (TouchSurface::Touch& touch)
        {
            auto& startPos = touchStartPositions.getValue (touch);

            if (touch.isTouchStart)
                startPos = { touch.x, touch.y };

            touch.startX = startPos.x;
            touch.startY = startPos.y;
        }

        void handlePacketACK (BlocksProtocol::TopologyIndex deviceIndex,
                              BlocksProtocol::PacketCounter counter)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
            {
                detector.handleSharedDataACK (deviceID, counter);
                updateApiPing (deviceID);
            }
        }

        void handleFirmwareUpdateACK (BlocksProtocol::TopologyIndex deviceIndex,
                                      BlocksProtocol::FirmwareUpdateACKCode resultCode,
                                      BlocksProtocol::FirmwareUpdateACKDetail resultDetail)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
            {
                detector.handleFirmwareUpdateACK (deviceID, (uint8) resultCode.get(), (uint32) resultDetail.get());
                updateApiPing (deviceID);
            }
        }

        void handleConfigUpdateMessage (BlocksProtocol::TopologyIndex deviceIndex,
                                        int32 item, int32 value, int32 min, int32 max)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleConfigUpdateMessage (deviceID, item, value, min, max);
        }

        void handleConfigSetMessage (BlocksProtocol::TopologyIndex deviceIndex,
                                     int32 item, int32 value)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleConfigSetMessage (deviceID, item, value);
        }

        void handleConfigFactorySyncEndMessage (BlocksProtocol::TopologyIndex deviceIndex)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleConfigFactorySyncEndMessage (deviceID);
        }

        void handleConfigFactorySyncResetMessage (BlocksProtocol::TopologyIndex deviceIndex)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleConfigFactorySyncResetMessage (deviceID);
        }

        void handleLogMessage (BlocksProtocol::TopologyIndex deviceIndex, const String& message)
        {
            if (auto deviceID = getDeviceIDFromMessageIndex (deviceIndex))
                detector.handleLogMessage (deviceID, message);
        }

        //==============================================================================
        template <typename PacketBuilder>
        bool sendMessageToDevice (const PacketBuilder& builder) const
        {
            if (deviceConnection->sendMessageToDevice (builder.getData(), (size_t) builder.size()))
            {
               #if DUMP_BANDWIDTH_STATS
                registerBytesOut (builder.size());
               #endif
                return true;
            }

            return false;
        }

        DeviceConnection* getDeviceConnection()
        {
            return deviceConnection.get();
        }

        juce::Array<DeviceInfo> getCurrentDeviceInfo()
        {
            auto blocks = currentDeviceInfo;
            blocks.removeIf ([this] (DeviceInfo& info) { return ! isApiConnected (info.uid); });
            return blocks;
        }

        juce::Array<BlockDeviceConnection> getCurrentDeviceConnections()
        {
            auto connections = currentDeviceConnections;
            connections.removeIf ([this] (BlockDeviceConnection& c) { return ! isApiConnected (c.device1) || ! isApiConnected (c.device2); });
            return connections;
        }

        Detector& detector;
        juce::String deviceName;

        static constexpr double pingTimeoutSeconds = 6.0;

    private:
        //==============================================================================
        juce::Array<DeviceInfo> currentDeviceInfo;
        juce::Array<BlockDeviceConnection> currentDeviceConnections;
        std::unique_ptr<DeviceConnection> deviceConnection;

        juce::Array<BlocksProtocol::DeviceStatus> incomingTopologyDevices, currentTopologyDevices;
        juce::Array<BlocksProtocol::DeviceConnection> incomingTopologyConnections;

        juce::CriticalSection incomingPacketLock;
        juce::Array<juce::MemoryBlock> incomingPackets;

        struct TouchStart { float x, y; };
        TouchList<TouchStart> touchStartPositions;

        //==============================================================================
        juce::Time lastTopologyRequestTime, lastTopologyReceiveTime;
        int numTopologyRequestsSent = 0;

        void scheduleNewTopologyRequest()
        {
            numTopologyRequestsSent = 0;
            lastTopologyReceiveTime = juce::Time();
            lastTopologyRequestTime = juce::Time::getCurrentTime();
        }

        void sendTopologyRequest()
        {
            ++numTopologyRequestsSent;
            lastTopologyRequestTime = juce::Time::getCurrentTime();
            sendCommandMessage (0, BlocksProtocol::requestTopologyMessage);
        }

        void timerCallback() override
        {
            const auto now = juce::Time::getCurrentTime();

            if ((now > lastTopologyReceiveTime + juce::RelativeTime::seconds (30.0))
                && now > lastTopologyRequestTime + juce::RelativeTime::seconds (1.0)
                && numTopologyRequestsSent < 4)
                sendTopologyRequest();

            checkApiTimeouts (now);
            startApiModeOnConnectedBlocks();
        }

        bool failedToGetTopology() const noexcept
        {
            return numTopologyRequestsSent > 4 && lastTopologyReceiveTime == juce::Time();
        }

        bool sendCommandMessage (BlocksProtocol::TopologyIndex deviceIndex, uint32 commandID) const
        {
            BlocksProtocol::HostPacketBuilder<64> p;
            p.writePacketSysexHeaderBytes (deviceIndex);
            p.deviceControlMessage (commandID);
            p.writePacketSysexFooter();
            return sendMessageToDevice (p);
        }

        //==============================================================================
        struct BlockPingTime
        {
            Block::UID blockUID;
            juce::Time lastPing;
        };

        juce::Array<BlockPingTime> blockPings;

        void updateApiPing (Block::UID uid)
        {
            const auto now = juce::Time::getCurrentTime();

            if (auto* ping = getPing (uid))
            {
                LOG_PING ("Ping: " << uid << " " << now.formatted ("%Mm %Ss"));
                ping->lastPing = now;
            }
            else
            {
                LOG_CONNECTIVITY ("API Connected " << uid);
                blockPings.add ({ uid, now });
                detector.handleTopologyChange();
            }
        }

        BlockPingTime* getPing (Block::UID uid)
        {
            for (auto& ping : blockPings)
                if (uid == ping.blockUID)
                    return &ping;

            return nullptr;
        }

        void removeDeviceInfo (Block::UID uid)
        {
            currentDeviceInfo.removeIf ([uid] (DeviceInfo& info) { return uid == info.uid; });
        }

        bool isApiConnected (Block::UID uid)
        {
            return getPing (uid) != nullptr;
        }

        void forceApiDisconnected (Block::UID uid)
        {
            if (isApiConnected (uid))
            {
                // Clear all known API connections and broadcast an empty topology,
                // as DNA blocks connected to the restarting block may be offline.
                LOG_CONNECTIVITY ("API Disconnected " << uid << ", re-probing topology");
                currentDeviceInfo.clearQuick();
                blockPings.clearQuick();
                detector.handleTopologyChange();
                scheduleNewTopologyRequest();
            }
        }

        void checkApiTimeouts (juce::Time now)
        {
            const auto timedOut = [this, now] (BlockPingTime& ping)
            {
                if (ping.lastPing >= now - juce::RelativeTime::seconds (pingTimeoutSeconds))
                    return false;

                LOG_CONNECTIVITY ("Ping timeout: " << ping.blockUID);
                removeDeviceInfo (ping.blockUID);
                return true;
            };

            if (blockPings.removeIf (timedOut) > 0)
            {
                scheduleNewTopologyRequest();
                detector.handleTopologyChange();
            }
        }

        void startApiModeOnConnectedBlocks()
        {
            for (auto& info : currentDeviceInfo)
            {
                if (! isApiConnected (info.uid))
                {
                    LOG_CONNECTIVITY ("API Try " << info.uid);
                    sendCommandMessage (info.index, BlocksProtocol::endAPIMode);
                    sendCommandMessage (info.index, BlocksProtocol::beginAPIMode);
                }
            }
        }

        //==============================================================================
        Block::UID getDeviceIDFromIndex (BlocksProtocol::TopologyIndex index) const noexcept
        {
            for (auto& d : currentDeviceInfo)
                if (d.index == index)
                    return d.uid;

            return {};
        }

        Block::UID getDeviceIDFromMessageIndex (BlocksProtocol::TopologyIndex index) noexcept
        {
            const auto uid = getDeviceIDFromIndex (index);

            // re-request topology if we get an event from an unknown block
            if (uid == Block::UID())
                scheduleNewTopologyRequest();

            return uid;
        }

        juce::Array<BlockDeviceConnection> getArrayOfConnections (const juce::Array<BlocksProtocol::DeviceConnection>& connections)
        {
            juce::Array<BlockDeviceConnection> result;

            for (auto&& c : connections)
            {
                BlockDeviceConnection dc;
                dc.device1 = getDeviceIDFromIndex (c.device1);
                dc.device2 = getDeviceIDFromIndex (c.device2);

                if (dc.device1 <= 0 || dc.device2 <= 0)
                    continue;

                dc.connectionPortOnDevice1 = convertConnectionPort (dc.device1, c.port1);
                dc.connectionPortOnDevice2 = convertConnectionPort (dc.device2, c.port2);

                result.add (dc);
            }

            return result;
        }

        Block::ConnectionPort convertConnectionPort (Block::UID uid, BlocksProtocol::ConnectorPort p) noexcept
        {
            if (auto* info = getDeviceInfoFromUID (uid))
                return BlocksProtocol::BlockDataSheet (info->serial).convertPortIndexToConnectorPort (p);

            jassertfalse;
            return { Block::ConnectionPort::DeviceEdge::north, 0 };
        }

        //==============================================================================
        void handleIncomingMessage (const void* data, size_t dataSize)
        {
            juce::MemoryBlock mb (data, dataSize);

            {
                const juce::ScopedLock sl (incomingPacketLock);
                incomingPackets.add (std::move (mb));
            }

            triggerAsyncUpdate();

           #if DUMP_BANDWIDTH_STATS
            registerBytesIn ((int) dataSize);
           #endif
        }

        void handleAsyncUpdate() override
        {
            juce::Array<juce::MemoryBlock> packets;
            packets.ensureStorageAllocated (32);

            {
                const juce::ScopedLock sl (incomingPacketLock);
                incomingPackets.swapWith (packets);
            }

            for (auto& packet : packets)
            {
                auto data = static_cast<const uint8*> (packet.getData());

                BlocksProtocol::HostPacketDecoder<ConnectedDeviceGroup>
                    ::processNextPacket (*this, *data, data + 1, (int) packet.getSize() - 1);
            }
        }

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ConnectedDeviceGroup)
    };

    //==============================================================================
    /** This is the main singleton object that keeps track of connected blocks */
    struct Detector   : public juce::ReferenceCountedObject,
                        private juce::Timer
    {
        Detector()  : defaultDetector (new MIDIDeviceDetector()), deviceDetector (*defaultDetector)
        {
            startTimer (10);
        }

        Detector (DeviceDetector& dd)  : deviceDetector (dd)
        {
            startTimer (10);
        }

        ~Detector()
        {
            jassert (activeTopologySources.isEmpty());
        }

        using Ptr = juce::ReferenceCountedObjectPtr<Detector>;

        static Detector::Ptr getDefaultDetector()
        {
            auto& d = getDefaultDetectorPointer();

            if (d == nullptr)
                d = new Detector();

            return d;
        }

        static Detector::Ptr& getDefaultDetectorPointer()
        {
            static Detector::Ptr defaultDetector;
            return defaultDetector;
        }

        void detach (PhysicalTopologySource* pts)
        {
            activeTopologySources.removeAllInstancesOf (pts);

            if (activeTopologySources.isEmpty())
            {
                for (auto& b : currentTopology.blocks)
                    if (auto bi = BlockImplementation::getFrom (*b))
                        bi->sendCommandMessage (BlocksProtocol::endAPIMode);

                currentTopology = {};
                lastTopology = {};

                auto& d = getDefaultDetectorPointer();

                if (d != nullptr && d->getReferenceCount() == 2)
                    getDefaultDetectorPointer() = nullptr;
            }
        }

        bool isConnected (Block::UID deviceID) const noexcept
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED // This method must only be called from the message thread!

            for (auto&& b : currentTopology.blocks)
                if (b->uid == deviceID)
                    return true;

            return false;
        }

        const BlocksProtocol::DeviceStatus* getLastStatus (Block::UID deviceID) const noexcept
        {
            for (auto d : connectedDeviceGroups)
                if (auto status = d->getLastStatus (deviceID))
                    return status;

            return nullptr;
        }

        void handleTopologyChange()
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

            {
                juce::Array<DeviceInfo> newDeviceInfo;
                juce::Array<BlockDeviceConnection> newDeviceConnections;

                for (auto d : connectedDeviceGroups)
                {
                    newDeviceInfo.addArray (d->getCurrentDeviceInfo());
                    newDeviceConnections.addArray (d->getCurrentDeviceConnections());
                }

                for (int i = currentTopology.blocks.size(); --i >= 0;)
                {
                    auto currentBlock = currentTopology.blocks.getUnchecked (i);

                    auto newDeviceIter = std::find_if (newDeviceInfo.begin(), newDeviceInfo.end(),
                                                       [&] (DeviceInfo& info) { return info.uid == currentBlock->uid; });

                    auto* blockImpl = BlockImplementation::getFrom (*currentBlock);

                    if (newDeviceIter == newDeviceInfo.end())
                    {
                        if (blockImpl != nullptr)
                            blockImpl->markDisconnected();

                        disconnectedBlocks.addIfNotAlreadyThere (currentTopology.blocks.removeAndReturn (i).get());
                    }
                    else
                    {
                        if (blockImpl != nullptr && blockImpl->wasPowerCycled())
                        {
                            blockImpl->resetPowerCycleFlag();
                            blockImpl->markReconnected (newDeviceIter->version, newDeviceIter->name, newDeviceIter->isMaster);
                        }

                        updateCurrentBlockInfo (currentBlock, *newDeviceIter);
                    }
                }

                static const int maxBlocksToSave = 100;

                if (disconnectedBlocks.size() > maxBlocksToSave)
                    disconnectedBlocks.removeRange (0, 2 * (disconnectedBlocks.size() - maxBlocksToSave));

                for (auto& info : newDeviceInfo)
                    if (info.serial.isValid() && ! containsBlockWithUID (currentTopology.blocks, getBlockUIDFromSerialNumber (info.serial)))
                        addBlock (info);

                currentTopology.connections.swapWith (newDeviceConnections);
            }

            broadcastTopology();
        }

        void notifyBlockIsRestarting (Block::UID deviceID)
        {
            for (auto& group : connectedDeviceGroups)
                group->notifyBlockIsRestarting (deviceID);
        }

        void handleSharedDataACK (Block::UID deviceID, uint32 packetCounter) const
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

            if (auto* bi = getBlockImplementationWithUID (deviceID))
                bi->handleSharedDataACK (packetCounter);
        }

        void handleFirmwareUpdateACK (Block::UID deviceID, uint8 resultCode, uint32 resultDetail)
        {
            if (auto* bi = getBlockImplementationWithUID (deviceID))
                bi->handleFirmwareUpdateACK (resultCode, resultDetail);
        }

        void handleConfigUpdateMessage (Block::UID deviceID, int32 item, int32 value, int32 min, int32 max)
        {
            if (auto* bi = getBlockImplementationWithUID (deviceID))
                bi->handleConfigUpdateMessage (item, value, min, max);
        }

        void notifyBlockOfConfigChange (BlockImplementation& bi, uint32 item)
        {
            if (auto configChangedCallback = bi.configChangedCallback)
            {
                if (item >= bi.getMaxConfigIndex())
                    configChangedCallback (bi, {}, item);
                else
                    configChangedCallback (bi, bi.getLocalConfigMetaData (item), item);
            }
        }

        void handleConfigSetMessage (Block::UID deviceID, int32 item, int32 value)
        {
            if (auto* bi = getBlockImplementationWithUID (deviceID))
            {
                bi->handleConfigSetMessage (item, value);
                notifyBlockOfConfigChange (*bi, uint32 (item));
            }
        }

        void handleConfigFactorySyncEndMessage (Block::UID deviceID)
        {
            if (auto* bi = getBlockImplementationWithUID (deviceID))
                notifyBlockOfConfigChange (*bi, bi->getMaxConfigIndex());
        }

        void handleConfigFactorySyncResetMessage (Block::UID deviceID)
        {
            if (auto* bi = getBlockImplementationWithUID (deviceID))
                bi->resetConfigListActiveStatus();
        }

        void handleLogMessage (Block::UID deviceID, const String& message) const
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

            if (auto* bi = getBlockImplementationWithUID (deviceID))
                bi->handleLogMessage (message);
        }

        void handleButtonChange (Block::UID deviceID, Block::Timestamp timestamp, uint32 buttonIndex, bool isDown) const
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

            if (auto* bi = getBlockImplementationWithUID (deviceID))
            {
                bi->pingFromDevice();

                if (isPositiveAndBelow (buttonIndex, bi->getButtons().size()))
                    if (auto* cbi = dynamic_cast<ControlButtonImplementation*> (bi->getButtons().getUnchecked (int (buttonIndex))))
                        cbi->broadcastButtonChange (timestamp, bi->modelData.buttons[(int) buttonIndex].type, isDown);
            }
        }

        void handleTouchChange (Block::UID deviceID, const TouchSurface::Touch& touchEvent)
        {
            JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

            auto block = currentTopology.getBlockWithUID (deviceID);
            if (block != nullptr)
            {
                if (auto* surface = dynamic_cast<TouchSurfaceImplementation*> (block->getTouchSurface()))
                {
                    TouchSurface::Touch scaledEvent (touchEvent);

                    scaledEvent.x      *= block->getWidth();
                    scaledEvent.y      *= block->getHeight();
                    scaledEvent.startX *= block->getWidth();
                    scaledEvent.startY *= block->getHeight();

                    surface->broadcastTouchChange (scaledEvent);
                }
            }
        }

        void cancelAllActiveTouches() noexcept
        {
            for (auto& block : currentTopology.blocks)
                if (auto* surface = block->getTouchSurface())
                    surface->cancelAllActiveTouches();
        }

        void handleCustomMessage (Block::UID deviceID, Block::Timestamp timestamp, const int32* data)
        {
            if (auto* bi = getBlockImplementationWithUID (deviceID))
                bi->handleCustomMessage (timestamp, data);
        }

        //==============================================================================
        int getIndexFromDeviceID (Block::UID deviceID) const noexcept
        {
            for (auto* c : connectedDeviceGroups)
            {
                auto index = c->getIndexFromDeviceID (deviceID);

                if (index >= 0)
                    return index;
            }

            return -1;
        }

        template <typename PacketBuilder>
        bool sendMessageToDevice (Block::UID deviceID, const PacketBuilder& builder) const
        {
            for (auto* c : connectedDeviceGroups)
                if (c->getIndexFromDeviceID (deviceID) >= 0)
                    return c->sendMessageToDevice (builder);

            return false;
        }

        static Detector* getFrom (Block& b) noexcept
        {
            if (auto* bi = BlockImplementation::getFrom (b))
                return (bi->detector);

            jassertfalse;
            return nullptr;
        }

        DeviceConnection* getDeviceConnectionFor (const Block& b)
        {
            for (const auto& d : connectedDeviceGroups)
            {
                for (const auto& info : d->getCurrentDeviceInfo())
                {
                    if (info.uid == b.uid)
                        return d->getDeviceConnection();
                }
            }

            return nullptr;
        }

        const DeviceConnection* getDeviceConnectionFor (const Block& b) const
        {
            for (const auto& d : connectedDeviceGroups)
            {
                for (const auto& info : d->getCurrentDeviceInfo())
                {
                    if (info.uid == b.uid)
                        return d->getDeviceConnection();
                }
            }

            return nullptr;
        }

        std::unique_ptr<MIDIDeviceDetector> defaultDetector;
        DeviceDetector& deviceDetector;

        juce::Array<PhysicalTopologySource*> activeTopologySources;

        BlockTopology currentTopology, lastTopology;
        juce::ReferenceCountedArray<Block, CriticalSection> disconnectedBlocks;

    private:
        void timerCallback() override
        {
            startTimer (1500);

            auto detectedDevices = deviceDetector.scanForDevices();

            handleDevicesRemoved (detectedDevices);
            handleDevicesAdded (detectedDevices);
        }

        void handleDevicesRemoved (const juce::StringArray& detectedDevices)
        {
            bool anyDevicesRemoved = false;

            for (int i = connectedDeviceGroups.size(); --i >= 0;)
            {
                if (! connectedDeviceGroups.getUnchecked(i)->isStillConnected (detectedDevices))
                {
                    connectedDeviceGroups.remove (i);
                    anyDevicesRemoved = true;
                }
            }

            if (anyDevicesRemoved)
                handleTopologyChange();
        }

        void handleDevicesAdded (const juce::StringArray& detectedDevices)
        {
            for (const auto& devName : detectedDevices)
            {
                if (! hasDeviceFor (devName))
                {
                    if (auto d = deviceDetector.openDevice (detectedDevices.indexOf (devName)))
                    {
                        connectedDeviceGroups.add (new ConnectedDeviceGroup (*this, devName, d));
                    }
                }
            }
        }

        bool hasDeviceFor (const juce::String& devName) const
        {
            for (auto d : connectedDeviceGroups)
                if (d->deviceName == devName)
                    return true;

            return false;
        }

        void addBlock (DeviceInfo info)
        {
            if (! reactivateBlockIfKnown (info))
                addNewBlock (info);
        }

        bool reactivateBlockIfKnown (DeviceInfo info)
        {
            const auto uid = getBlockUIDFromSerialNumber (info.serial);

            for (int i = disconnectedBlocks.size(); --i >= 0;)
            {
                if (uid != disconnectedBlocks.getUnchecked (i)->uid)
                    continue;

                auto block = disconnectedBlocks.removeAndReturn (i);

                if (auto* blockImpl = BlockImplementation::getFrom (*block))
                {
                    blockImpl->markReconnected (info.version, info.name, info.isMaster);
                    currentTopology.blocks.add (block);
                    return true;
                }
            }

            return false;
        }

        void addNewBlock (DeviceInfo info)
        {
            currentTopology.blocks.add (new BlockImplementation (info.serial, *this, info.version,
                                                                 info.name, info.isMaster));
        }

        void updateCurrentBlockInfo (Block::Ptr blockToUpdate, DeviceInfo& updatedInfo)
        {
            if (versionNumberChanged (updatedInfo, blockToUpdate->versionNumber))
                setVersionNumberForBlock (updatedInfo, *blockToUpdate);

            if (nameIsValid (updatedInfo))
                setNameForBlock (updatedInfo, *blockToUpdate);

            if (updatedInfo.isMaster != blockToUpdate->isMasterBlock())
                BlockImplementation::getFrom (*blockToUpdate)->setToMaster (updatedInfo.isMaster);
        }

        BlockImplementation* getBlockImplementationWithUID (Block::UID deviceID) const noexcept
        {
            if (auto&& block = currentTopology.getBlockWithUID (deviceID))
                return BlockImplementation::getFrom (*block);

            return nullptr;
        }

        juce::OwnedArray<ConnectedDeviceGroup> connectedDeviceGroups;

        //==============================================================================
        void broadcastTopology()
        {
            if (currentTopology != lastTopology)
            {
                lastTopology = currentTopology;

                BlocksTraverser traverser;
                traverser.traverseBlockArray (currentTopology);

                for (auto* d : activeTopologySources)
                    d->listeners.call ([] (TopologySource::Listener& l) { l.topologyChanged(); });

               #if DUMP_TOPOLOGY
                dumpTopology (lastTopology);
               #endif
            }
        }

        JUCE_DECLARE_WEAK_REFERENCEABLE (Detector)
        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Detector)
    };

    //==============================================================================
    /** This is a friend of the BlocksImplementation that will scan and set the
        physical positions of the blocks */
    struct BlocksTraverser
    {
        void traverseBlockArray (const BlockTopology& topology)
        {
            juce::Array<Block::UID> visited;

            for (auto& block : topology.blocks)
            {
                if (block->isMasterBlock() && ! visited.contains (block->uid))
                {
                    if (auto* bi = dynamic_cast<BlockImplementation*> (block))
                    {
                        bi->masterUID = {};
                        bi->position = {};
                        bi->rotation = 0;
                    }

                    layoutNeighbours (*block, topology, block->uid, visited);
                }
            }
        }

        // returns the distance from corner clockwise
        int getUnitForIndex (Block::Ptr block, Block::ConnectionPort::DeviceEdge edge, int index)
        {
            if (block->getType() == Block::seaboardBlock)
            {
                if (edge == Block::ConnectionPort::DeviceEdge::north)
                {
                    if (index == 0) return 1;
                    if (index == 1) return 4;
                }
                else if (edge != Block::ConnectionPort::DeviceEdge::south)
                {
                    return 1;
                }
            }

            if (edge == Block::ConnectionPort::DeviceEdge::south)
                return block->getWidth() - (index + 1);

            if (edge == Block::ConnectionPort::DeviceEdge::west)
                return block->getHeight() - (index + 1);

            return index;
        }

        // returns how often north needs to rotate by 90 degrees
        int getRotationForEdge (Block::ConnectionPort::DeviceEdge edge)
        {
            switch (edge)
            {
                case Block::ConnectionPort::DeviceEdge::north:  return 0;
                case Block::ConnectionPort::DeviceEdge::east:   return 1;
                case Block::ConnectionPort::DeviceEdge::south:  return 2;
                case Block::ConnectionPort::DeviceEdge::west:   return 3;
            }

            jassertfalse;
            return 0;
        }

        void layoutNeighbours (Block::Ptr block, const BlockTopology& topology,
                               Block::UID masterUid, juce::Array<Block::UID>& visited)
        {
            visited.add (block->uid);

            for (auto& connection : topology.connections)
            {
                if ((connection.device1 == block->uid && ! visited.contains (connection.device2))
                     || (connection.device2 == block->uid && ! visited.contains (connection.device1)))
                {
                    const auto theirUid = connection.device1 == block->uid ? connection.device2 : connection.device1;
                    const auto neighbourPtr = topology.getBlockWithUID (theirUid);

                    if (auto* neighbour = dynamic_cast<BlockImplementation*> (neighbourPtr.get()))
                    {
                        const auto  myBounds    = block->getBlockAreaWithinLayout();
                        const auto& myPort      = connection.device1 == block->uid ? connection.connectionPortOnDevice1 : connection.connectionPortOnDevice2;
                        const auto& theirPort   = connection.device1 == block->uid ? connection.connectionPortOnDevice2 : connection.connectionPortOnDevice1;
                        const auto  myOffset    = getUnitForIndex (block, myPort.edge, myPort.index);
                        const auto  theirOffset = getUnitForIndex (neighbourPtr, theirPort.edge, theirPort.index);

                        neighbour->masterUID = masterUid;
                        neighbour->rotation = (2 + block->getRotation()
                                                 + getRotationForEdge (myPort.edge)
                                                 - getRotationForEdge (theirPort.edge)) % 4;

                        Point<int> delta;
                        const auto theirBounds = neighbour->getBlockAreaWithinLayout();

                        switch ((block->getRotation() + getRotationForEdge (myPort.edge)) % 4)
                        {
                            case 0: // over me
                                delta = { myOffset - (theirBounds.getWidth() - (theirOffset + 1)), -theirBounds.getHeight() };
                                break;
                            case 1: // right of me
                                delta = { myBounds.getWidth(), myOffset - (theirBounds.getHeight() - (theirOffset + 1)) };
                                break;
                            case 2: // under me
                                delta = { (myBounds.getWidth() - (myOffset + 1)) - theirOffset, myBounds.getHeight() };
                                break;
                            case 3: // left of me
                                delta = { -theirBounds.getWidth(), (myBounds.getHeight() - (myOffset + 1)) - theirOffset };
                                break;
                        }

                        neighbour->position = myBounds.getPosition() + delta;
                    }

                    layoutNeighbours (neighbourPtr, topology, masterUid, visited);
                }
            }
        }
    };

    //==============================================================================
    struct BlockImplementation  : public Block,
                                  private MIDIDeviceConnection::Listener,
                                  private Timer
    {
        BlockImplementation (const BlocksProtocol::BlockSerialNumber& serial,
                             Detector& detectorToUse,
                             BlocksProtocol::VersionNumber version,
                             BlocksProtocol::BlockName blockName,
                             bool isMasterBlock)
            : Block (juce::String ((const char*) serial.serial,   sizeof (serial.serial)),
                     juce::String ((const char*) version.version, version.length),
                     juce::String ((const char*) blockName.name,  blockName.length)),
              modelData (serial),
              remoteHeap (modelData.programAndHeapSize),
              detector (&detectorToUse),
              isMaster (isMasterBlock)
        {
            if (modelData.hasTouchSurface)
                touchSurface.reset (new TouchSurfaceImplementation (*this));

            int i = 0;

            for (auto&& b : modelData.buttons)
                controlButtons.add (new ControlButtonImplementation (*this, i++, b));

            if (modelData.lightGridWidth > 0 && modelData.lightGridHeight > 0)
                ledGrid.reset (new LEDGridImplementation (*this));

            for (auto&& s : modelData.statusLEDs)
                statusLights.add (new StatusLightImplementation (*this, s));

            updateMidiConnectionListener();
        }

        ~BlockImplementation()
        {
            if (listenerToMidiConnection != nullptr)
            {
                config.setDeviceComms (nullptr);
                listenerToMidiConnection->removeListener (this);
            }
        }

        void markDisconnected()
        {
            if (auto surface = dynamic_cast<TouchSurfaceImplementation*> (touchSurface.get()))
                surface->disableTouchSurface();
        }

        void markReconnected (BlocksProtocol::VersionNumber newVersion, BlocksProtocol::BlockName newName, bool master)
        {
            versionNumber = getVersionString (newVersion);
            name = getNameString (newName);
            isMaster = master;

            setProgram (nullptr);
            remoteHeap.resetDeviceStateToUnknown();

            if (auto surface = dynamic_cast<TouchSurfaceImplementation*> (touchSurface.get()))
                surface->activateTouchSurface();

            updateMidiConnectionListener();
        }

        void setToMaster (bool shouldBeMaster)
        {
            isMaster = shouldBeMaster;
        }

        void updateMidiConnectionListener()
        {
            if (detector == nullptr)
                return;

            listenerToMidiConnection = dynamic_cast<MIDIDeviceConnection*> (detector->getDeviceConnectionFor (*this));

            if (listenerToMidiConnection != nullptr)
                listenerToMidiConnection->addListener (this);

            config.setDeviceComms (listenerToMidiConnection);
        }

        Type getType() const override                                   { return modelData.apiType; }
        juce::String getDeviceDescription() const override              { return modelData.description; }
        int getWidth() const override                                   { return modelData.widthUnits; }
        int getHeight() const override                                  { return modelData.heightUnits; }
        float getMillimetersPerUnit() const override                    { return 47.0f; }
        bool isHardwareBlock() const override                           { return true; }
        juce::Array<Block::ConnectionPort> getPorts() const override    { return modelData.ports; }
        bool isConnected() const override                               { return detector && detector->isConnected (uid); }
        bool isMasterBlock() const override                             { return isMaster; }
        Block::UID getConnectedMasterUID() const override               { return masterUID; }
        int getRotation() const override                                { return rotation; }

        Rectangle<int> getBlockAreaWithinLayout() const override
        {
            if (rotation % 2 == 0)
                return { position.getX(), position.getY(), modelData.widthUnits, modelData.heightUnits };

            return { position.getX(), position.getY(), modelData.heightUnits, modelData.widthUnits };
        }

        TouchSurface* getTouchSurface() const override                  { return touchSurface.get(); }
        LEDGrid* getLEDGrid() const override                            { return ledGrid.get(); }

        LEDRow* getLEDRow() override
        {
            if (ledRow == nullptr && modelData.numLEDRowLEDs > 0)
                ledRow.reset (new LEDRowImplementation (*this));

            return ledRow.get();
        }

        juce::Array<ControlButton*> getButtons() const override
        {
            juce::Array<ControlButton*> result;
            result.addArray (controlButtons);
            return result;
        }

        juce::Array<StatusLight*> getStatusLights() const override
        {
            juce::Array<StatusLight*> result;
            result.addArray (statusLights);
            return result;
        }

        float getBatteryLevel() const override
        {
            if (detector == nullptr)
                return 0.0f;

            if (auto status = detector->getLastStatus (uid))
                return status->batteryLevel.toUnipolarFloat();

            return 0.0f;
        }

        bool isBatteryCharging() const override
        {
            if (detector == nullptr)
                return false;

            if (auto status = detector->getLastStatus (uid))
                return status->batteryCharging.get() != 0;

            return false;
        }

        bool supportsGraphics() const override
        {
            return false;
        }

        int getDeviceIndex() const noexcept
        {
            if (detector == nullptr)
                return -1;

            return isConnected() ? detector->getIndexFromDeviceID (uid) : -1;
        }

        template <typename PacketBuilder>
        bool sendMessageToDevice (const PacketBuilder& builder)
        {
            if (detector != nullptr)
            {
                lastMessageSendTime = juce::Time::getCurrentTime();
                return detector->sendMessageToDevice (uid, builder);
            }

            return false;
        }

        bool sendCommandMessage (uint32 commandID)
        {
            return buildAndSendPacket<64> ([commandID] (BlocksProtocol::HostPacketBuilder<64>& p)
                                           { return p.deviceControlMessage (commandID); });
        }

        void handleCustomMessage (Block::Timestamp, const int32* data)
        {
            ProgramEventMessage m;

            for (uint32 i = 0; i < BlocksProtocol::numProgramMessageInts; ++i)
                m.values[i] = data[i];

            programEventListeners.call ([&] (ProgramEventListener& l) { l.handleProgramEvent (*this, m); });
        }

        static BlockImplementation* getFrom (Block& b) noexcept
        {
            jassert (dynamic_cast<BlockImplementation*> (&b) != nullptr);
            return dynamic_cast<BlockImplementation*> (&b);
        }

        //==============================================================================
        std::function<void(const String&)> logger;

        void setLogger (std::function<void(const String&)> newLogger) override
        {
            logger = newLogger;
        }

        void handleLogMessage (const String& message) const
        {
            if (logger != nullptr)
                logger (message);
        }

        //==============================================================================
        juce::Result setProgram (Program* newProgram) override
        {
            if (newProgram != nullptr && program.get() == newProgram)
            {
                jassertfalse;
                return juce::Result::ok();
            }

            stopTimer();

            {
                std::unique_ptr<Program> p (newProgram);

                if (program != nullptr
                    && newProgram != nullptr
                    && program->getLittleFootProgram() == newProgram->getLittleFootProgram())
                    return juce::Result::ok();

                std::swap (program, p);
            }

            programSize = 0;
            isProgramLoaded = shouldSaveProgramAsDefault = false;

            if (program == nullptr)
            {
                remoteHeap.clear();
                return juce::Result::ok();
            }

            littlefoot::Compiler compiler;
            compiler.addNativeFunctions (PhysicalTopologySource::getStandardLittleFootFunctions());

            const auto err = compiler.compile (program->getLittleFootProgram(), 512, program->getSearchPaths());

            if (err.failed())
                return err;

            DBG ("Compiled littlefoot program, space needed: "
                 << (int) compiler.getCompiledProgram().getTotalSpaceNeeded() << " bytes");

            if (compiler.getCompiledProgram().getTotalSpaceNeeded() > getMemorySize())
                return Result::fail ("Program too large!");

            const auto size = (size_t) compiler.compiledObjectCode.size();
            programSize = (uint32) size;

            remoteHeap.resetDataRangeToUnknown (0, remoteHeap.blockSize);
            remoteHeap.clear();
            remoteHeap.sendChanges (*this, true);

            remoteHeap.resetDataRangeToUnknown (0, (uint32) size);
            remoteHeap.setBytes (0, compiler.compiledObjectCode.begin(), size);
            remoteHeap.sendChanges (*this, true);

            this->resetConfigListActiveStatus();

            if (auto changeCallback = this->configChangedCallback)
                changeCallback (*this, {}, this->getMaxConfigIndex());

            startTimer (20);

            return juce::Result::ok();
        }

        Program* getProgram() const override                                        { return program.get(); }

        void sendProgramEvent (const ProgramEventMessage& message) override
        {
            static_assert (sizeof (ProgramEventMessage::values) == 4 * BlocksProtocol::numProgramMessageInts,
                           "Need to keep the internal and external messages structures the same");

            if (remoteHeap.isProgramLoaded())
            {
                buildAndSendPacket<128> ([&message] (BlocksProtocol::HostPacketBuilder<128>& p)
                                         { return p.addProgramEventMessage (message.values); });
            }
        }

        void timerCallback() override
        {
            if (remoteHeap.isFullySynced() && remoteHeap.isProgramLoaded())
            {
                isProgramLoaded = true;
                stopTimer();

                if (shouldSaveProgramAsDefault)
                    doSaveProgramAsDefault();

                if (programLoadedCallback != nullptr)
                    programLoadedCallback (*this);
            }
            else
            {
                startTimer (100);
            }
        }

        void saveProgramAsDefault() override
        {
            shouldSaveProgramAsDefault = true;

            if (! isTimerRunning() && isProgramLoaded)
                doSaveProgramAsDefault();
        }

        uint32 getMemorySize() override
        {
            return modelData.programAndHeapSize;
        }

        uint32 getHeapMemorySize() override
        {
            jassert (isPositiveAndNotGreaterThan (programSize, modelData.programAndHeapSize));
            return modelData.programAndHeapSize - programSize;
        }

        void setDataByte (size_t offset, uint8 value) override
        {
            remoteHeap.setByte (programSize + offset, value);
        }

        void setDataBytes (size_t offset, const void* newData, size_t num) override
        {
            remoteHeap.setBytes (programSize + offset, static_cast<const uint8*> (newData), num);
        }

        void setDataBits (uint32 startBit, uint32 numBits, uint32 value) override
        {
            remoteHeap.setBits (programSize * 8 + startBit, numBits, value);
        }

        uint8 getDataByte (size_t offset) override
        {
            return remoteHeap.getByte (programSize + offset);
        }

        void handleSharedDataACK (uint32 packetCounter) noexcept
        {
            pingFromDevice();
            remoteHeap.handleACKFromDevice (*this, packetCounter);
        }

        bool sendFirmwareUpdatePacket (const uint8* data, uint8 size, std::function<void (uint8, uint32)> callback) override
        {
            firmwarePacketAckCallback = {};

            if (buildAndSendPacket<256> ([data, size] (BlocksProtocol::HostPacketBuilder<256>& p)
                                         { return p.addFirmwareUpdatePacket (data, size); }))
            {
                firmwarePacketAckCallback = callback;
                return true;
            }

            return false;
        }

        void handleFirmwareUpdateACK (uint8 resultCode, uint32 resultDetail)
        {
            if (firmwarePacketAckCallback != nullptr)
            {
                firmwarePacketAckCallback (resultCode, resultDetail);
                firmwarePacketAckCallback = {};
            }
        }

        void handleConfigUpdateMessage (int32 item, int32 value, int32 min, int32 max)
        {
            config.handleConfigUpdateMessage (item, value, min, max);
        }

        void handleConfigSetMessage(int32 item, int32 value)
        {
            config.handleConfigSetMessage (item, value);
        }

        void pingFromDevice()
        {
            lastMessageReceiveTime = juce::Time::getCurrentTime();
        }

        void addDataInputPortListener (DataInputPortListener* listener) override
        {
            Block::addDataInputPortListener (listener);

            if (auto midiInput = getMidiInput())
                midiInput->start();
        }

        void sendMessage (const void* message, size_t messageSize) override
        {
            if (auto midiOutput = getMidiOutput())
                midiOutput->sendMessageNow ({ message, (int) messageSize });
        }

        void handleTimerTick()
        {
            if (ledGrid != nullptr)
                if (auto renderer = ledGrid->getRenderer())
                    renderer->renderLEDGrid (*ledGrid);

            remoteHeap.sendChanges (*this, false);

            if (lastMessageSendTime < juce::Time::getCurrentTime() - juce::RelativeTime::milliseconds (pingIntervalMs))
                sendCommandMessage (BlocksProtocol::ping);
        }

        //==============================================================================
        int32 getLocalConfigValue (uint32 item) override
        {
            initialiseDeviceIndexAndConnection();
            return config.getItemValue ((BlocksProtocol::ConfigItemId) item);
        }

        void setLocalConfigValue (uint32 item, int32 value) override
        {
            initialiseDeviceIndexAndConnection();
            config.setItemValue ((BlocksProtocol::ConfigItemId) item, value);
        }

        void setLocalConfigRange (uint32 item, int32 min, int32 max) override
        {
            initialiseDeviceIndexAndConnection();
            config.setItemMin ((BlocksProtocol::ConfigItemId) item, min);
            config.setItemMax ((BlocksProtocol::ConfigItemId) item, max);
        }

        void setLocalConfigItemActive (uint32 item, bool isActive) override
        {
            initialiseDeviceIndexAndConnection();
            config.setItemActive ((BlocksProtocol::ConfigItemId) item, isActive);
        }

        bool isLocalConfigItemActive (uint32 item) override
        {
            initialiseDeviceIndexAndConnection();
            return config.getItemActive ((BlocksProtocol::ConfigItemId) item);
        }

        uint32 getMaxConfigIndex() override
        {
            return uint32 (BlocksProtocol::maxConfigIndex);
        }

        bool isValidUserConfigIndex (uint32 item) override
        {
            return item >= (uint32) BlocksProtocol::ConfigItemId::user0
                && item < (uint32) (BlocksProtocol::ConfigItemId::user0 + numberOfUserConfigs);
        }

        ConfigMetaData getLocalConfigMetaData (uint32 item) override
        {
            initialiseDeviceIndexAndConnection();
            return config.getMetaData ((BlocksProtocol::ConfigItemId) item);
        }

        void requestFactoryConfigSync() override
        {
            initialiseDeviceIndexAndConnection();
            config.requestFactoryConfigSync();
        }

        void resetConfigListActiveStatus() override
        {
            config.resetConfigListActiveStatus();
        }

        void setConfigChangedCallback (std::function<void(Block&, const ConfigMetaData&, uint32)> configChanged) override
        {
            configChangedCallback = std::move (configChanged);
        }

        void setProgramLoadedCallback (std::function<void(Block&)> programLoaded) override
        {
            programLoadedCallback = std::move (programLoaded);
        }

        bool setName (const juce::String& newName) override
        {
            return buildAndSendPacket<128> ([&newName] (BlocksProtocol::HostPacketBuilder<128>& p)
                                            { return p.addSetBlockName (newName); });
        }

        void factoryReset() override
        {
            buildAndSendPacket<32> ([] (BlocksProtocol::HostPacketBuilder<32>& p)
                                    { return p.addFactoryReset(); });
        }

        void blockReset() override
        {
            if (buildAndSendPacket<32> ([] (BlocksProtocol::HostPacketBuilder<32>& p)
                                        { return p.addBlockReset(); }))
            {
                hasBeenPowerCycled = true;

                if (detector != nullptr)
                    detector->notifyBlockIsRestarting (uid);
            }
        }

        bool wasPowerCycled() const { return hasBeenPowerCycled; }
        void resetPowerCycleFlag()  { hasBeenPowerCycled = false; }

        //==============================================================================
        std::unique_ptr<TouchSurface> touchSurface;
        juce::OwnedArray<ControlButton> controlButtons;
        std::unique_ptr<LEDGridImplementation> ledGrid;
        std::unique_ptr<LEDRowImplementation> ledRow;
        juce::OwnedArray<StatusLight> statusLights;

        BlocksProtocol::BlockDataSheet modelData;

        MIDIDeviceConnection* listenerToMidiConnection = nullptr;

        static constexpr int pingIntervalMs = 400;

        static constexpr uint32 maxBlockSize = BlocksProtocol::padBlockProgramAndHeapSize;
        static constexpr uint32 maxPacketCounter = BlocksProtocol::PacketCounter::maxValue;
        static constexpr uint32 maxPacketSize = 200;

        using PacketBuilder = BlocksProtocol::HostPacketBuilder<maxPacketSize>;

        using RemoteHeapType = littlefoot::LittleFootRemoteHeap<BlockImplementation>;
        RemoteHeapType remoteHeap;

        WeakReference<Detector> detector;
        juce::Time lastMessageSendTime, lastMessageReceiveTime;

        BlockConfigManager config;
        std::function<void(Block&, const ConfigMetaData&, uint32)> configChangedCallback;

        std::function<void(Block&)> programLoadedCallback;

    private:
        std::unique_ptr<Program> program;
        uint32 programSize = 0;

        std::function<void(uint8, uint32)> firmwarePacketAckCallback;

        bool isMaster = false;
        Block::UID masterUID = {};

        Point<int> position;
        int rotation = 0;
        friend BlocksTraverser;

        bool isProgramLoaded = false;
        bool shouldSaveProgramAsDefault = false;
        bool hasBeenPowerCycled = false;

        void initialiseDeviceIndexAndConnection()
        {
            config.setDeviceIndex ((TopologyIndex) getDeviceIndex());
            config.setDeviceComms (listenerToMidiConnection);
        }

        const juce::MidiInput* getMidiInput() const
        {
            if (detector != nullptr)
                if (auto c = dynamic_cast<const MIDIDeviceConnection*> (detector->getDeviceConnectionFor (*this)))
                    return c->midiInput.get();

            jassertfalse;
            return nullptr;
        }

        juce::MidiInput* getMidiInput()
        {
            return const_cast<juce::MidiInput*> (static_cast<const BlockImplementation&>(*this).getMidiInput());
        }

        const juce::MidiOutput* getMidiOutput() const
        {
            if (detector != nullptr)
                if (auto c = dynamic_cast<const MIDIDeviceConnection*> (detector->getDeviceConnectionFor (*this)))
                    return c->midiOutput.get();

            jassertfalse;
            return nullptr;
        }

        juce::MidiOutput* getMidiOutput()
        {
            return const_cast<juce::MidiOutput*> (static_cast<const BlockImplementation&>(*this).getMidiOutput());
        }

        void handleIncomingMidiMessage (const juce::MidiMessage& message) override
        {
            dataInputPortListeners.call ([&] (DataInputPortListener& l) { l.handleIncomingDataPortMessage (*this, message.getRawData(),
                                                                                                           (size_t) message.getRawDataSize()); });
        }

        void connectionBeingDeleted (const MIDIDeviceConnection& c) override
        {
            jassert (listenerToMidiConnection == &c);
            juce::ignoreUnused (c);
            listenerToMidiConnection->removeListener (this);
            listenerToMidiConnection = nullptr;
            config.setDeviceComms (nullptr);
        }

        void doSaveProgramAsDefault()
        {
            sendCommandMessage (BlocksProtocol::saveProgramAsDefault);
        }

        template<int packetBytes, typename PacketBuilderFn>
        bool buildAndSendPacket (PacketBuilderFn buildFn)
        {
            auto index = getDeviceIndex();

            if (index < 0)
            {
                jassertfalse;
                return false;
            }

            BlocksProtocol::HostPacketBuilder<packetBytes> p;
            p.writePacketSysexHeaderBytes ((BlocksProtocol::TopologyIndex) index);

            if (! buildFn (p))
                return false;

            p.writePacketSysexFooter();
            return sendMessageToDevice (p);
        }

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (BlockImplementation)
    };

    //==============================================================================
    struct LEDRowImplementation  : public LEDRow,
                                   private Timer
    {
        LEDRowImplementation (BlockImplementation& b) : LEDRow (b)
        {
            startTimer (300);
        }

        void setButtonColour (uint32 index, LEDColour colour)
        {
            if (index < 10)
            {
                colours[index] = colour;
                flush();
            }
        }

        int getNumLEDs() const override
        {
            return static_cast<const BlockImplementation&> (block).modelData.numLEDRowLEDs;
        }

        void setLEDColour (int index, LEDColour colour) override
        {
            if ((uint32) index < 15u)
            {
                colours[10 + index] = colour;
                flush();
            }
        }

        void setOverlayColour (LEDColour colour) override
        {
            colours[25] = colour;
            flush();
        }

        void resetOverlayColour() override
        {
            setOverlayColour ({});
        }

    private:
        LEDColour colours[26];

        void timerCallback() override
        {
            stopTimer();
            loadProgramOntoBlock();
            flush();
        }

        void loadProgramOntoBlock()
        {
            if (block.getProgram() == nullptr)
            {
                auto err = block.setProgram (new DefaultLEDGridProgram (block));

                if (err.failed())
                {
                    DBG (err.getErrorMessage());
                    jassertfalse;
                }
            }
        }

        void flush()
        {
            if (block.getProgram() != nullptr)
                for (uint32 i = 0; i < (uint32) numElementsInArray (colours); ++i)
                    write565Colour (16 * i, colours[i]);
        }

        void write565Colour (uint32 bitIndex, LEDColour colour)
        {
            block.setDataBits (bitIndex,      5, colour.getRed()   >> 3);
            block.setDataBits (bitIndex + 5,  6, colour.getGreen() >> 2);
            block.setDataBits (bitIndex + 11, 5, colour.getBlue()  >> 3);
        }

        struct DefaultLEDGridProgram  : public Block::Program
        {
            DefaultLEDGridProgram (Block& b) : Block::Program (b) {}

            juce::String getLittleFootProgram() override
            {
                /*  Data format:

                    0:  10 x 5-6-5 bits for button LED RGBs
                    20: 15 x 5-6-5 bits for LED row colours
                    50:  1 x 5-6-5 bits for LED row overlay colour
                */
                return R"littlefoot(

                #heapsize: 128

                int getColour (int bitIndex)
                {
                    return makeARGB (255,
                                     getHeapBits (bitIndex,      5) << 3,
                                     getHeapBits (bitIndex + 5,  6) << 2,
                                     getHeapBits (bitIndex + 11, 5) << 3);
                }

                int getButtonColour (int index)
                {
                    return getColour (16 * index);
                }

                int getLEDColour (int index)
                {
                    if (getHeapInt (50))
                        return getColour (50 * 8);

                    return getColour (20 * 8 + 16 * index);
                }

                void repaint()
                {
                    for (int x = 0; x < 15; ++x)
                        fillPixel (getLEDColour (x), x, 0);

                    for (int i = 0; i < 10; ++i)
                        fillPixel (getButtonColour (i), i, 1);
                }

                void handleMessage (int p1, int p2) {}

                )littlefoot";
            }
        };

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LEDRowImplementation)
    };

    //==============================================================================
    struct TouchSurfaceImplementation  : public TouchSurface,
                                         private juce::Timer
    {
        TouchSurfaceImplementation (BlockImplementation& b)  : TouchSurface (b), blockImpl (b)
        {
            activateTouchSurface();
        }

        ~TouchSurfaceImplementation()
        {
            disableTouchSurface();
        }

        void activateTouchSurface()
        {
            startTimer (500);
        }

        void disableTouchSurface()
        {
            stopTimer();
        }

        int getNumberOfKeywaves() const noexcept override
        {
            return blockImpl.modelData.numKeywaves;
        }

        void broadcastTouchChange (const TouchSurface::Touch& touchEvent)
        {
            auto& status = touches.getValue (touchEvent);

            // Fake a touch end if we receive a duplicate touch-start with no preceding touch-end (ie: comms error)
            if (touchEvent.isTouchStart && status.isActive)
                killTouch (touchEvent, status, juce::Time::getMillisecondCounter());

            // Fake a touch start if we receive an unexpected event with no matching start event. (ie: comms error)
            if (! touchEvent.isTouchStart && ! status.isActive)
            {
                TouchSurface::Touch t (touchEvent);
                t.isTouchStart = true;
                t.isTouchEnd = false;

                if (t.zVelocity <= 0)  t.zVelocity = status.lastStrikePressure;
                if (t.zVelocity <= 0)  t.zVelocity = t.z;
                if (t.zVelocity <= 0)  t.zVelocity = 0.9f;

                listeners.call ([&] (TouchSurface::Listener& l) { l.touchChanged (*this, t); });
            }

            // Normal handling:
            status.lastEventTime = juce::Time::getMillisecondCounter();
            status.isActive = ! touchEvent.isTouchEnd;

            if (touchEvent.isTouchStart)
                status.lastStrikePressure = touchEvent.zVelocity;

            listeners.call ([&] (TouchSurface::Listener& l) { l.touchChanged (*this, touchEvent); });
        }

        void timerCallback() override
        {
            // Find touches that seem to have become stuck, and fake a touch-end for them..
            static const uint32 touchTimeOutMs = 500;

            for (auto& t : touches)
            {
                auto& status = t.value;
                auto now = juce::Time::getMillisecondCounter();

                if (status.isActive && now > status.lastEventTime + touchTimeOutMs)
                    killTouch (t.touch, status, now);
            }
        }

        struct TouchStatus
        {
            uint32 lastEventTime = 0;
            float lastStrikePressure = 0;
            bool isActive = false;
        };

        void killTouch (const TouchSurface::Touch& touch, TouchStatus& status, uint32 timeStamp) noexcept
        {
            jassert (status.isActive);

            TouchSurface::Touch killTouch (touch);

            killTouch.z                 = 0;
            killTouch.xVelocity         = 0;
            killTouch.yVelocity         = 0;
            killTouch.zVelocity         = -1.0f;
            killTouch.eventTimestamp    = timeStamp;
            killTouch.isTouchStart      = false;
            killTouch.isTouchEnd        = true;

            listeners.call ([&] (TouchSurface::Listener& l) { l.touchChanged (*this, killTouch); });

            status.isActive = false;
        }

        void cancelAllActiveTouches() noexcept override
        {
            const auto now = juce::Time::getMillisecondCounter();

            for (auto& t : touches)
                if (t.value.isActive)
                    killTouch (t.touch, t.value, now);

            touches.clear();
        }

        BlockImplementation& blockImpl;
        TouchList<TouchStatus> touches;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (TouchSurfaceImplementation)
    };

    //==============================================================================
    struct ControlButtonImplementation  : public ControlButton
    {
        ControlButtonImplementation (BlockImplementation& b, int index, BlocksProtocol::BlockDataSheet::ButtonInfo info)
            : ControlButton (b), blockImpl (b), buttonInfo (info), buttonIndex (index)
        {
        }

        ~ControlButtonImplementation()
        {
        }

        ButtonFunction getType() const override         { return buttonInfo.type; }
        juce::String getName() const override           { return BlocksProtocol::getButtonNameForFunction (buttonInfo.type); }
        float getPositionX() const override             { return buttonInfo.x; }
        float getPositionY() const override             { return buttonInfo.y; }

        bool hasLight() const override                  { return blockImpl.isControlBlock(); }

        bool setLightColour (LEDColour colour) override
        {
            if (hasLight())
            {
                if (auto row = blockImpl.ledRow.get())
                {
                    row->setButtonColour ((uint32) buttonIndex, colour);
                    return true;
                }
            }

            return false;
        }

        void broadcastButtonChange (Block::Timestamp timestamp, ControlButton::ButtonFunction button, bool isDown)
        {
            if (button == buttonInfo.type)
            {
                if (wasDown == isDown)
                    sendButtonChangeToListeners (timestamp, ! isDown);

                sendButtonChangeToListeners (timestamp, isDown);
                wasDown = isDown;
            }
        }

        void sendButtonChangeToListeners (Block::Timestamp timestamp, bool isDown)
        {
            if (isDown)
                listeners.call ([&] (ControlButton::Listener& l) { l.buttonPressed (*this, timestamp); });
            else
                listeners.call ([&] (ControlButton::Listener& l) { l.buttonReleased (*this, timestamp); });
        }

        BlockImplementation& blockImpl;
        BlocksProtocol::BlockDataSheet::ButtonInfo buttonInfo;
        int buttonIndex;
        bool wasDown = false;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ControlButtonImplementation)
    };


    //==============================================================================
    struct StatusLightImplementation  : public StatusLight
    {
        StatusLightImplementation (Block& b, BlocksProtocol::BlockDataSheet::StatusLEDInfo i)  : StatusLight (b), info (i)
        {
        }

        juce::String getName() const override               { return info.name; }

        bool setColour (LEDColour newColour) override
        {
            // XXX TODO!
            juce::ignoreUnused (newColour);
            return false;
        }

        BlocksProtocol::BlockDataSheet::StatusLEDInfo info;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (StatusLightImplementation)
    };

    //==============================================================================
    struct LEDGridImplementation  : public LEDGrid
    {
        LEDGridImplementation (BlockImplementation& b)  : LEDGrid (b), blockImpl (b)
        {
        }

        int getNumColumns() const override      { return blockImpl.modelData.lightGridWidth; }
        int getNumRows() const override         { return blockImpl.modelData.lightGridHeight; }

        BlockImplementation& blockImpl;

        JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LEDGridImplementation)
    };

    //==============================================================================
   #if DUMP_TOPOLOGY
    static juce::String idToSerialNum (const BlockTopology& topology, Block::UID uid)
    {
        for (auto* b : topology.blocks)
            if (b->uid == uid)
                return b->serialNumber;

        return "???";
    }

    static juce::String portEdgeToString (Block::ConnectionPort port)
    {
        switch (port.edge)
        {
            case Block::ConnectionPort::DeviceEdge::north: return "north";
            case Block::ConnectionPort::DeviceEdge::south: return "south";
            case Block::ConnectionPort::DeviceEdge::east:  return "east";
            case Block::ConnectionPort::DeviceEdge::west:  return "west";
        }

        return {};
    }

    static juce::String portToString (Block::ConnectionPort port)
    {
        return portEdgeToString (port) + "_" + juce::String (port.index);
    }

    static void dumpTopology (const BlockTopology& topology)
    {
        MemoryOutputStream m;

        m << "=============================================================================" << newLine
          << "Topology:  " << topology.blocks.size() << " device(s)" << newLine
          << newLine;

        int index = 0;

        for (auto block : topology.blocks)
        {
            m << "Device " << index++ << (block->isMasterBlock() ? ":  (MASTER)" : ":") << newLine;

            m << "  Description: " << block->getDeviceDescription() << newLine
              << "  Serial: " << block->serialNumber << newLine;

            if (auto bi = BlockImplementation::getFrom (*block))
                m << "  Short address: " << (int) bi->getDeviceIndex() << newLine;

            m << "  Battery level: " + juce::String (juce::roundToInt (100.0f * block->getBatteryLevel())) + "%" << newLine
              << "  Battery charging: " + juce::String (block->isBatteryCharging() ? "y" : "n") << newLine
              << "  Width: " << block->getWidth() << newLine
              << "  Height: " << block->getHeight() << newLine
              << "  Millimeters per unit: " << block->getMillimetersPerUnit() << newLine
              << newLine;
        }

        for (auto& connection : topology.connections)
        {
            m << idToSerialNum (topology, connection.device1)
              << ":" << portToString (connection.connectionPortOnDevice1)
              << "  <->  "
              << idToSerialNum (topology, connection.device2)
              << ":" << portToString (connection.connectionPortOnDevice2) << newLine;
        }

        m << "=============================================================================" << newLine;

        Logger::outputDebugString (m.toString());
    }
   #endif
};

//==============================================================================
struct PhysicalTopologySource::DetectorHolder  : private juce::Timer
{
    DetectorHolder (PhysicalTopologySource& pts)
        : topologySource (pts),
          detector (Internal::Detector::getDefaultDetector())
    {
        startTimerHz (30);
    }

    DetectorHolder (PhysicalTopologySource& pts, DeviceDetector& dd)
        : topologySource (pts),
          detector (new Internal::Detector (dd))
    {
        startTimerHz (30);
    }

    void timerCallback() override
    {
        if (! topologySource.hasOwnServiceTimer())
            handleTimerTick();
    }

    void handleTimerTick()
    {
        for (auto& b : detector->currentTopology.blocks)
            if (auto bi = Internal::BlockImplementation::getFrom (*b))
                bi->handleTimerTick();
    }

    PhysicalTopologySource& topologySource;
    Internal::Detector::Ptr detector;
};

//==============================================================================
PhysicalTopologySource::PhysicalTopologySource (bool startDetached)
{
    if (! startDetached)
        setActive (true);
}

PhysicalTopologySource::PhysicalTopologySource (DeviceDetector& detectorToUse, bool startDetached)
    : customDetector (&detectorToUse)
{
    if (! startDetached)
        setActive (true);
}

PhysicalTopologySource::~PhysicalTopologySource()
{
    setActive (false);
}

void PhysicalTopologySource::setActive (bool shouldBeActive)
{
    JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

    if (isActive() == shouldBeActive)
        return;

    if (shouldBeActive)
    {
        if (customDetector == nullptr)
            detector = std::make_unique<DetectorHolder>(*this);
        else
            detector = std::make_unique<DetectorHolder>(*this, *customDetector);

        detector->detector->activeTopologySources.add (this);
    }
    else
    {
        detector->detector->detach (this);
        detector.reset();
    }

    listeners.call ([](TopologySource::Listener& l){ l.topologyChanged(); });
}

bool PhysicalTopologySource::isActive() const
{
    return detector != nullptr;
}

bool PhysicalTopologySource::isLockedFromOutside() const
{
    if (detector != nullptr && detector->detector != nullptr)
        return detector->detector->deviceDetector.isLockedFromOutside();

    return false;
}

BlockTopology PhysicalTopologySource::getCurrentTopology() const
{
    JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED // This method must only be called from the message thread!

    if (detector != nullptr)
        return detector->detector->currentTopology;

    return {};
}

void PhysicalTopologySource::cancelAllActiveTouches() noexcept
{
    if (detector != nullptr)
        detector->detector->cancelAllActiveTouches();
}

bool PhysicalTopologySource::hasOwnServiceTimer() const     { return false; }
void PhysicalTopologySource::handleTimerTick()
{
    if (detector != nullptr)
        detector->handleTimerTick();
}

PhysicalTopologySource::DeviceConnection::DeviceConnection() {}
PhysicalTopologySource::DeviceConnection::~DeviceConnection() {}

PhysicalTopologySource::DeviceDetector::DeviceDetector() {}
PhysicalTopologySource::DeviceDetector::~DeviceDetector() {}

const char* const* PhysicalTopologySource::getStandardLittleFootFunctions() noexcept
{
    return BlocksProtocol::ledProgramLittleFootFunctions;
}

template <typename ListType>
static bool collectionsMatch (const ListType& list1, const ListType& list2) noexcept
{
    if (list1.size() != list2.size())
        return false;

    for (auto&& b : list1)
        if (! list2.contains (b))
            return false;

    return true;
}

bool BlockTopology::operator== (const BlockTopology& other) const noexcept
{
    return collectionsMatch (connections, other.connections) && collectionsMatch (blocks, other.blocks);
}

bool BlockTopology::operator!= (const BlockTopology& other) const noexcept
{
    return ! operator== (other);
}

bool BlockDeviceConnection::operator== (const BlockDeviceConnection& other) const noexcept
{
    return (device1 == other.device1 && device2 == other.device2
             && connectionPortOnDevice1 == other.connectionPortOnDevice1
             && connectionPortOnDevice2 == other.connectionPortOnDevice2)
        || (device1 == other.device2 && device2 == other.device1
             && connectionPortOnDevice1 == other.connectionPortOnDevice2
             && connectionPortOnDevice2 == other.connectionPortOnDevice1);
}

bool BlockDeviceConnection::operator!= (const BlockDeviceConnection& other) const noexcept
{
    return ! operator== (other);
}

} // namespace juce
