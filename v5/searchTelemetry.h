#pragma once

#ifndef ASSEMBLY_ENABLE_TELEMETRY

inline constexpr bool searchTelemetryCompiled = false;
inline bool searchTelemetryEnabled = false;

#else

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <type_traits>

#ifdef __linux__
    #include <cerrno>
    #include <cstdlib>
    #include <cstring>
    #include <fcntl.h>
    #include <unistd.h>
#endif

enum class SearchTelemetryPhase : size_t
{
    inputSetup,
    initialEnumeration,
    dagConversion,
    assemblySearch,
    output,
    count
};

struct SearchTelemetryCounters
{
    uint64_t retainedMaskAttempts = 0;
    uint64_t retainedMasks = 0;
    uint64_t duplicateMaskAttempts = 0;
    uint64_t rejectedMasks = 0;
    uint64_t matchingVisits = 0;

    uint64_t canonicalisationCalls = 0;
    uint64_t canonicalisationMaskCacheHits = 0;
    uint64_t canonicalisationMaskCacheMisses = 0;
    uint64_t canonicalClassInsertions = 0;
    uint64_t canonicalClassReuses = 0;
    uint64_t vf2Calls = 0;
    uint64_t vf2Matches = 0;

    uint64_t residualDecompositionRequests = 0;
    uint64_t residualCacheEligibleRequests = 0;
    uint64_t residualCacheSmallMoleculeBypasses = 0;
    uint64_t residualCacheWideMoleculeBypasses = 0;
    uint64_t residualCacheSmallResidualBypasses = 0;
    uint64_t residualCacheFirstOccurrenceBypasses = 0;
    uint64_t residualCacheRuntimeDisabledBypasses = 0;
    uint64_t residualCacheLookups = 0;
    uint64_t residualCacheHits = 0;
    uint64_t residualCacheMisses = 0;
    uint64_t residualCacheAdmissions = 0;

    uint64_t assemblyCacheLookups = 0;
    uint64_t assemblyCacheHits = 0;
    uint64_t assemblyCacheMisses = 0;
    uint64_t assemblyCachePrunedHits = 0;
    uint64_t assemblyCacheUpdatedHits = 0;

    uint64_t pairBoundCacheLookups = 0;
    uint64_t pairBoundCacheHits = 0;
    uint64_t pairBoundCacheMisses = 0;
};

struct ProcessMemorySnapshot
{
    bool residentAvailable = false;
    bool virtualAvailable = false;
    uint64_t residentKiB = 0;
    uint64_t residentHighWaterKiB = 0;
    uint64_t virtualKiB = 0;
    uint64_t virtualPeakKiB = 0;
};

struct SearchTelemetryPhaseStats
{
    uint64_t clockTicks = 0;
    uint64_t startResidentKiB = 0;
    uint64_t peakResidentKiB = 0;
    uint64_t endResidentKiB = 0;
    uint64_t startVirtualKiB = 0;
    uint64_t endVirtualKiB = 0;
    size_t activations = 0;
    bool startResidentAvailable = false;
    bool endResidentAvailable = false;
    bool startVirtualAvailable = false;
    bool endVirtualAvailable = false;
    bool exactResidentPeak = false;
};

struct SearchTelemetryState
{
    SearchTelemetryCounters counters;
    std::array<
        SearchTelemetryPhaseStats,
        static_cast<size_t>(SearchTelemetryPhase::count)
    > phases;
    SearchTelemetryPhase currentPhase = SearchTelemetryPhase::count;
    clock_t phaseStart = 0;
    size_t processedAtoms = 0;
    size_t processedEdges = 0;
    size_t activeMaskWords = 0;
    bool residualCacheEligible = false;
    bool active = false;
    ProcessMemorySnapshot finalMemory;
};

inline constexpr bool searchTelemetryCompiled = true;

inline bool searchTelemetryEnabled = false;
inline ASSEMBLYCPP_SEARCH_LOCAL SearchTelemetryState searchTelemetry;

inline const char* searchTelemetryPhaseName(SearchTelemetryPhase phase)
{
    switch (phase)
    {
        case SearchTelemetryPhase::inputSetup: return "input_setup";
        case SearchTelemetryPhase::initialEnumeration:
            return "initial_enumeration";
        case SearchTelemetryPhase::dagConversion: return "dag_conversion";
        case SearchTelemetryPhase::assemblySearch: return "assembly_search";
        case SearchTelemetryPhase::output: return "output";
        case SearchTelemetryPhase::count: break;
    }
    return "unknown";
}

inline uint64_t telemetryClockDifference(clock_t start, clock_t end)
{
    const clock_t clockError = static_cast<clock_t>(-1);
    if (start == clockError || end == clockError) return 0;
    using UnsignedClock = std::make_unsigned_t<clock_t>;
    return static_cast<uint64_t>(
        static_cast<UnsignedClock>(end) - static_cast<UnsignedClock>(start)
    );
}

#ifdef __linux__
inline bool parseProcStatusValue(
    const char *buffer,
    const char *end,
    const char *label,
    uint64_t &value
)
{
    const size_t labelLength = std::strlen(label);
    const char *position = buffer;
    while (position < end)
    {
        const char *lineEnd = static_cast<const char *>(
            std::memchr(position, '\n', static_cast<size_t>(end - position))
        );
        if (lineEnd == nullptr) lineEnd = end;
        if (
            static_cast<size_t>(lineEnd - position) >= labelLength &&
            std::memcmp(position, label, labelLength) == 0
        )
        {
            const char *number = position + labelLength;
            while (number < lineEnd && (*number == ' ' || *number == '\t'))
                ++number;
            errno = 0;
            char *parsedEnd = nullptr;
            const unsigned long long parsed = std::strtoull(
                number,
                &parsedEnd,
                10
            );
            if (
                errno != 0 ||
                parsedEnd == number ||
                parsedEnd > lineEnd
            ) return false;
            value = static_cast<uint64_t>(parsed);
            return true;
        }
        position = lineEnd < end ? lineEnd + 1 : end;
    }
    return false;
}

inline ProcessMemorySnapshot readProcessMemorySnapshot()
{
    ProcessMemorySnapshot result;
    const int descriptor = open("/proc/self/status", O_RDONLY | O_CLOEXEC);
    if (descriptor < 0) return result;

    std::array<char, 16384> buffer{};
    size_t used = 0;
    while (used + 1 < buffer.size())
    {
        const ssize_t count = read(
            descriptor,
            buffer.data() + used,
            buffer.size() - used - 1
        );
        if (count == 0) break;
        if (count < 0)
        {
            if (errno == EINTR) continue;
            close(descriptor);
            return result;
        }
        used += static_cast<size_t>(count);
    }
    close(descriptor);
    const char *end = buffer.data() + used;
    result.residentAvailable =
        parseProcStatusValue(
            buffer.data(), end, "VmRSS:", result.residentKiB
        ) &&
        parseProcStatusValue(
            buffer.data(), end, "VmHWM:", result.residentHighWaterKiB
        );
    result.virtualAvailable =
        parseProcStatusValue(
            buffer.data(), end, "VmSize:", result.virtualKiB
        ) &&
        parseProcStatusValue(
            buffer.data(), end, "VmPeak:", result.virtualPeakKiB
        );
    return result;
}

inline bool resetProcessResidentHighWaterMark()
{
    const int descriptor = open("/proc/self/clear_refs", O_WRONLY | O_CLOEXEC);
    if (descriptor < 0) return false;
    constexpr char resetCommand[] = "5\n";
    ssize_t written = 0;
    do
    {
        written = write(descriptor, resetCommand, sizeof(resetCommand) - 1);
    } while (written < 0 && errno == EINTR);
    const bool closed = close(descriptor) == 0;
    return written == static_cast<ssize_t>(sizeof(resetCommand) - 1) && closed;
}
#else
inline ProcessMemorySnapshot readProcessMemorySnapshot()
{
    return {};
}

inline bool resetProcessResidentHighWaterMark()
{
    return false;
}
#endif

inline void finishSearchTelemetryPhase()
{
    if (
        !searchTelemetry.active ||
        searchTelemetry.currentPhase == SearchTelemetryPhase::count
    ) return;

    const clock_t now = clock();
    SearchTelemetryPhaseStats &stats = searchTelemetry.phases[
        static_cast<size_t>(searchTelemetry.currentPhase)
    ];
    stats.clockTicks += telemetryClockDifference(
        searchTelemetry.phaseStart,
        now
    );
    ProcessMemorySnapshot memory = readProcessMemorySnapshot();
    stats.endResidentAvailable = memory.residentAvailable;
    if (memory.residentAvailable)
    {
        stats.endResidentKiB = memory.residentKiB;
        if (stats.exactResidentPeak)
        {
            if (memory.residentHighWaterKiB > stats.peakResidentKiB)
                stats.peakResidentKiB = memory.residentHighWaterKiB;
        }
    }
    else stats.exactResidentPeak = false;
    stats.endVirtualAvailable = memory.virtualAvailable;
    if (memory.virtualAvailable)
    {
        stats.endVirtualKiB = memory.virtualKiB;
    }
}

inline void setSearchTelemetryPhase(SearchTelemetryPhase phase)
{
    if (
        !searchTelemetryCompiled ||
        !searchTelemetryEnabled ||
        !searchTelemetry.active ||
        phase == SearchTelemetryPhase::count
    ) return;
    if (searchTelemetry.currentPhase == phase) return;
    finishSearchTelemetryPhase();

    const bool exactResidentPeak = resetProcessResidentHighWaterMark();
    ProcessMemorySnapshot memory = readProcessMemorySnapshot();
    SearchTelemetryPhaseStats &stats = searchTelemetry.phases[
        static_cast<size_t>(phase)
    ];
    const bool firstActivation = stats.activations == 0;
    ++stats.activations;
    stats.exactResidentPeak =
        (firstActivation || stats.exactResidentPeak) &&
        exactResidentPeak &&
        memory.residentAvailable;
    stats.endResidentAvailable = memory.residentAvailable;
    if (memory.residentAvailable)
    {
        if (firstActivation)
        {
            stats.startResidentAvailable = true;
            stats.startResidentKiB = memory.residentKiB;
        }
        stats.endResidentKiB = memory.residentKiB;
        if (exactResidentPeak)
        {
            if (memory.residentHighWaterKiB > stats.peakResidentKiB)
                stats.peakResidentKiB = memory.residentHighWaterKiB;
        }
    }
    stats.endVirtualAvailable = memory.virtualAvailable;
    if (memory.virtualAvailable)
    {
        if (firstActivation)
        {
            stats.startVirtualAvailable = true;
            stats.startVirtualKiB = memory.virtualKiB;
        }
        stats.endVirtualKiB = memory.virtualKiB;
    }
    searchTelemetry.currentPhase = phase;
    searchTelemetry.phaseStart = clock();
}

inline void resetSearchTelemetry()
{
    if (!searchTelemetryCompiled || !searchTelemetryEnabled) return;
    searchTelemetry = SearchTelemetryState{};
    searchTelemetry.active = true;
    setSearchTelemetryPhase(SearchTelemetryPhase::inputSetup);
}

inline void configureSearchTelemetryGraph(
    size_t atoms,
    size_t edges,
    size_t activeMaskWords,
    bool residualCacheEligible
)
{
    if (!searchTelemetryCompiled || !searchTelemetryEnabled) return;
    searchTelemetry.processedAtoms = atoms;
    searchTelemetry.processedEdges = edges;
    searchTelemetry.activeMaskWords = activeMaskWords;
    searchTelemetry.residualCacheEligible = residualCacheEligible;
}

inline void finaliseSearchTelemetry()
{
    if (
        !searchTelemetryCompiled ||
        !searchTelemetryEnabled ||
        !searchTelemetry.active
    ) return;
    finishSearchTelemetryPhase();
    searchTelemetry.finalMemory = readProcessMemorySnapshot();
    searchTelemetry.currentPhase = SearchTelemetryPhase::count;
    searchTelemetry.active = false;
}

inline void writeJsonRatio(
    std::ostream &output,
    uint64_t numerator,
    uint64_t denominator
)
{
    if (denominator == 0)
    {
        output << "null";
        return;
    }
    output << static_cast<double>(numerator) /
        static_cast<double>(denominator);
}

inline void writeJsonMemoryValue(
    std::ostream &output,
    bool available,
    uint64_t value
)
{
    if (available) output << value;
    else output << "null";
}

inline bool writeSearchTelemetry(const std::string &filename)
{
    if (!searchTelemetryCompiled || !searchTelemetryEnabled) return true;
    finaliseSearchTelemetry();

    std::ofstream output(filename);
    if (!output.is_open())
    {
        std::cerr << "error: could not open output file '"
                  << filename << "'\n";
        return false;
    }
    const SearchTelemetryCounters &counters = searchTelemetry.counters;
    const uint64_t canonicalLookups =
        counters.canonicalisationMaskCacheHits +
        counters.canonicalisationMaskCacheMisses;
    const uint64_t canonicalClassLookups =
        counters.canonicalClassInsertions + counters.canonicalClassReuses;
    const uint64_t residualLookupCount =
        counters.residualCacheHits + counters.residualCacheMisses;
    const uint64_t assemblyLookupCount =
        counters.assemblyCacheHits + counters.assemblyCacheMisses;
    const uint64_t pairBoundLookupCount =
        counters.pairBoundCacheHits + counters.pairBoundCacheMisses;

    bool anyActivatedPhase = false;
    bool allActivatedPhasePeaksExact = true;
    uint64_t overallResidentPeak = 0;
    for (const SearchTelemetryPhaseStats &phase : searchTelemetry.phases)
    {
        if (phase.activations == 0) continue;
        anyActivatedPhase = true;
        if (!phase.exactResidentPeak)
        {
            allActivatedPhasePeaksExact = false;
            continue;
        }
        overallResidentPeak = std::max(
            overallResidentPeak,
            phase.peakResidentKiB
        );
    }
    const bool overallResidentPeakAvailable =
        anyActivatedPhase && allActivatedPhasePeaksExact;

    output << "{\n"
           << "  \"schema_version\": 1,\n"
           << "  \"processed_graph\": {\n"
           << "    \"atoms\": " << searchTelemetry.processedAtoms << ",\n"
           << "    \"edges\": " << searchTelemetry.processedEdges << ",\n"
           << "    \"active_mask_words\": "
           << searchTelemetry.activeMaskWords << "\n"
           << "  },\n"
           << "  \"counters\": {\n"
           << "    \"retained_mask_attempts\": "
           << counters.retainedMaskAttempts << ",\n"
           << "    \"retained_masks\": " << counters.retainedMasks << ",\n"
           << "    \"duplicate_mask_attempts\": "
           << counters.duplicateMaskAttempts << ",\n"
           << "    \"rejected_masks\": " << counters.rejectedMasks << ",\n"
           << "    \"matching_visits\": " << counters.matchingVisits << ",\n"
           << "    \"canonicalisation_calls\": "
           << counters.canonicalisationCalls << ",\n"
           << "    \"vf2_calls\": " << counters.vf2Calls << ",\n"
           << "    \"vf2_matches\": " << counters.vf2Matches << "\n"
           << "  },\n"
           << "  \"caches\": {\n"
           << "    \"canonical_mask\": {\n"
           << "      \"hits\": "
           << counters.canonicalisationMaskCacheHits << ",\n"
           << "      \"misses\": "
           << counters.canonicalisationMaskCacheMisses << ",\n"
           << "      \"hit_rate\": ";
    writeJsonRatio(
        output,
        counters.canonicalisationMaskCacheHits,
        canonicalLookups
    );
    output << "\n    },\n"
           << "    \"canonical_class\": {\n"
           << "      \"insertions\": "
           << counters.canonicalClassInsertions << ",\n"
           << "      \"reuses\": " << counters.canonicalClassReuses << ",\n"
           << "      \"reuse_rate\": ";
    writeJsonRatio(output, counters.canonicalClassReuses, canonicalClassLookups);
    output << "\n    },\n"
           << "    \"residual_decomposition\": {\n"
           << "      \"eligible_for_processed_graph\": "
           << (searchTelemetry.residualCacheEligible ? "true" : "false")
           << ",\n"
           << "      \"requests\": "
           << counters.residualDecompositionRequests << ",\n"
           << "      \"eligible_requests\": "
           << counters.residualCacheEligibleRequests << ",\n"
           << "      \"small_molecule_bypasses\": "
           << counters.residualCacheSmallMoleculeBypasses << ",\n"
           << "      \"wide_molecule_bypasses\": "
           << counters.residualCacheWideMoleculeBypasses << ",\n"
           << "      \"small_residual_bypasses\": "
           << counters.residualCacheSmallResidualBypasses << ",\n"
           << "      \"first_occurrence_bypasses\": "
           << counters.residualCacheFirstOccurrenceBypasses << ",\n"
           << "      \"runtime_disabled_bypasses\": "
           << counters.residualCacheRuntimeDisabledBypasses << ",\n"
           << "      \"lookups\": " << counters.residualCacheLookups << ",\n"
           << "      \"hits\": " << counters.residualCacheHits << ",\n"
           << "      \"misses\": " << counters.residualCacheMisses << ",\n"
           << "      \"admissions\": "
           << counters.residualCacheAdmissions << ",\n"
           << "      \"lookup_hit_rate\": ";
    writeJsonRatio(output, counters.residualCacheHits, residualLookupCount);
    output << ",\n      \"request_hit_rate\": ";
    writeJsonRatio(
        output,
        counters.residualCacheHits,
        counters.residualDecompositionRequests
    );
    output << "\n    },\n"
           << "    \"assembly_state\": {\n"
           << "      \"lookups\": " << counters.assemblyCacheLookups << ",\n"
           << "      \"hits\": " << counters.assemblyCacheHits << ",\n"
           << "      \"misses\": " << counters.assemblyCacheMisses << ",\n"
           << "      \"pruned_hits\": "
           << counters.assemblyCachePrunedHits << ",\n"
           << "      \"updated_hits\": "
           << counters.assemblyCacheUpdatedHits << ",\n"
           << "      \"hit_rate\": ";
    writeJsonRatio(output, counters.assemblyCacheHits, assemblyLookupCount);
    output << "\n    },\n"
           << "    \"pair_bound\": {\n"
           << "      \"lookups\": " << counters.pairBoundCacheLookups << ",\n"
           << "      \"hits\": " << counters.pairBoundCacheHits << ",\n"
           << "      \"misses\": " << counters.pairBoundCacheMisses << ",\n"
           << "      \"hit_rate\": ";
    writeJsonRatio(output, counters.pairBoundCacheHits, pairBoundLookupCount);
    output << "\n    }\n"
           << "  },\n"
           << "  \"memory\": {\n"
           << "    \"method\": \"";
#ifdef __linux__
    output << "linux_proc_vmhwm_reset";
#else
    output << "unavailable";
#endif
    output << "\",\n"
           << "    \"phase_peaks_are_absolute_not_additive\": true,\n"
           << "    \"phase_peaks_complete\": "
           << (overallResidentPeakAvailable ? "true" : "false") << ",\n"
           << "    \"overall_peak_resident_kib\": ";
    writeJsonMemoryValue(
        output,
        overallResidentPeakAvailable,
        overallResidentPeak
    );
    output << ",\n    \"process_peak_virtual_kib\": ";
    writeJsonMemoryValue(
        output,
        searchTelemetry.finalMemory.virtualAvailable,
        searchTelemetry.finalMemory.virtualPeakKiB
    );
    output << ",\n    \"phases\": {\n";
    for (size_t i = 0; i < static_cast<size_t>(SearchTelemetryPhase::count); i++)
    {
        const SearchTelemetryPhase phase = static_cast<SearchTelemetryPhase>(i);
        const SearchTelemetryPhaseStats &stats = searchTelemetry.phases[i];
        const bool exact = stats.exactResidentPeak;
        output << "      \"" << searchTelemetryPhaseName(phase) << "\": {\n"
               << "        \"clock_ticks\": " << stats.clockTicks << ",\n"
               << "        \"activations\": " << stats.activations << ",\n"
               << "        \"start_rss_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.startResidentAvailable,
            stats.startResidentKiB
        );
        output << ",\n        \"peak_rss_kib\": ";
        writeJsonMemoryValue(output, exact, stats.peakResidentKiB);
        output << ",\n        \"end_rss_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.endResidentAvailable,
            stats.endResidentKiB
        );
        output << ",\n        \"start_virtual_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.startVirtualAvailable,
            stats.startVirtualKiB
        );
        output << ",\n        \"end_virtual_kib\": ";
        writeJsonMemoryValue(
            output,
            stats.endVirtualAvailable,
            stats.endVirtualKiB
        );
        output << "\n      }";
        if (i + 1 < static_cast<size_t>(SearchTelemetryPhase::count))
            output << ',';
        output << '\n';
    }
    output << "    }\n"
           << "  }\n"
           << "}\n";
    output.close();
    if (!output)
    {
        std::cerr << "error: could not write output file '"
                  << filename << "'\n";
        return false;
    }
    return true;
}

#endif
