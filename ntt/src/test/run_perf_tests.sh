#!/bin/bash

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Construct PROJECT_ROOT more robustly
# To reach 'Parallel-Programing-Project' from '.../ntt/src/test', we need to go up three levels.
PROJECT_ROOT="$(cd "${TEST_DIR}/../../../" && pwd)" # Go up three dirs and get absolute path

NTT_DATA_DIR_NAME=".nttdata"
NTT_DATA_DIR="${PROJECT_ROOT}/${NTT_DATA_DIR_NAME}"

CPP_TEST_HARNESS="${TEST_DIR}/test_ntt_openmp.cpp"
SERIAL_EXE="${TEST_DIR}/test_ntt_serial"
PARALLEL_EXE="${TEST_DIR}/test_ntt_parallel"

INCLUDE_DIR="${PROJECT_ROOT}/ntt/src/include" # For clarity in compiler command

ERROR_LOG_FILE="/tmp/ntt_test_error.log"

NUM_REPETITIONS=5 # Number of times to repeat each test for averaging

# Clean up previous error log if it exists
rm -f "${ERROR_LOG_FILE}"

# Debugging: Print resolved paths
echo "DEBUG: TEST_DIR is '${TEST_DIR}'"
echo "DEBUG: PROJECT_ROOT is '${PROJECT_ROOT}'"
echo "DEBUG: NTT_DATA_DIR_NAME is '${NTT_DATA_DIR_NAME}'"
echo "DEBUG: Attempting to check directory: '${NTT_DATA_DIR}'"

# Check if .nttdata directory exists
if [ ! -d "${NTT_DATA_DIR}" ]; then
    echo "Error: Test data directory not found at '${NTT_DATA_DIR}' (resolved from PROJECT_ROOT='${PROJECT_ROOT}')"
    echo "Please ensure the directory named '${NTT_DATA_DIR_NAME}' exists in the project root: '${PROJECT_ROOT}'"
    
    # Add more specific checks
    if [ ! -d "${PROJECT_ROOT}" ]; then
        echo "Critical Error: PROJECT_ROOT directory itself '${PROJECT_ROOT}' does not exist or is not a directory."
    else
        echo "Listing contents of PROJECT_ROOT ('${PROJECT_ROOT}') to help debug:"
        ls -la "${PROJECT_ROOT}"
    fi
    exit 1
fi


echo "--- Compiling C++ Test Harness (optimized for timing output) ---"
# Compile Serial Version
# The include path for basicMultiThread.h is relative from the .cpp file itself,
# so -I is not strictly needed if the #include is correct (e.g., "../include/thread/basicMultiThread.h")
# However, explicitly adding include path for ntt/src/include can be more robust if header structure changes slightly.
g++ -std=c++17 -O2 -I"${INCLUDE_DIR}" -o "${SERIAL_EXE}" "${CPP_TEST_HARNESS}"
if [ $? -ne 0 ]; then
    echo "Serial compilation failed!"
    exit 1
fi
echo "Serial version compiled: ${SERIAL_EXE}"

# Compile Parallel (OpenMP) Version
g++ -std=c++17 -O2 -fopenmp -I"${INCLUDE_DIR}" -o "${PARALLEL_EXE}" "${CPP_TEST_HARNESS}"
if [ $? -ne 0 ]; then
    echo "Parallel compilation failed!"
    exit 1
fi
echo "Parallel version compiled: ${PARALLEL_EXE}"
echo "------------------------------------"
echo ""

CPU_CORES=$(nproc --all 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) # Get core count, default to 4
THREAD_COUNTS=(1)
if [ "$CPU_CORES" -ge 2 ]; then THREAD_COUNTS+=(2); fi
if [ "$CPU_CORES" -ge 4 ]; then THREAD_COUNTS+=(4); fi
if (( CPU_CORES > 4 )); then # Add more threads if more cores
    if [ "$CPU_CORES" -ge 8 ]; then THREAD_COUNTS+=(8); fi
    if [ "$CPU_CORES" -ge 12 ]; then THREAD_COUNTS+=(12); fi # Example
    if [ "$CPU_CORES" -ge 16 ]; then THREAD_COUNTS+=(16); fi # Example
fi
# Ensure a test with actual core count if not already listed and not too large
if [[ ! " ${THREAD_COUNTS[@]} " =~ " ${CPU_CORES} " ]] && [ "$CPU_CORES" -le 16 ]; then
    THREAD_COUNTS+=($CPU_CORES)
    # Sort and unique to keep it tidy if CPU_CORES was already added by other conditions
    IFS=$'\n' THREAD_COUNTS=($(printf "%s\n" "${THREAD_COUNTS[@]}" | sort -un))
    unset IFS
fi


echo "--- Running Serial Tests (averaged over ${NUM_REPETITIONS} repetitions) ---"
OVERALL_SERIAL_STATUS="PASS"
# Test cases 0 through 3, skipping 4
TEST_CASES=(0 1 2 3)

for i in "${TEST_CASES[@]}"; do
    INPUT_FILE="${NTT_DATA_DIR}/${i}.in"
    OUTPUT_FILE="${NTT_DATA_DIR}/${i}.out"
    if [ ! -f "$INPUT_FILE" ] || [ ! -f "$OUTPUT_FILE" ]; then
        echo "Serial Test Case ${i} (${INPUT_FILE##*/}): Warning: Files not found. Skipping."
        OVERALL_SERIAL_STATUS="FAIL" # Or some other status like "SKIP"
        continue
    fi

    echo -n "Serial Test Case ${i} (${INPUT_FILE##*/}): "
    total_time_us=0
    run_status="PASS"
    min_time_us=-1
    max_time_us=0
    all_times_us=""

    for rep in $(seq 1 ${NUM_REPETITIONS}); do
        time_output=$("${SERIAL_EXE}" "${INPUT_FILE}" "${OUTPUT_FILE}" 2>>"${ERROR_LOG_FILE}")
        exit_code=$?
        if [ $exit_code -ne 0 ]; then
            run_status="FAIL (run ${rep}, code ${exit_code})"
            echo "" # Newline after prompt
            echo "Error during Serial Test Case ${i}, repetition ${rep}. Check ${ERROR_LOG_FILE}" >> /dev/stderr
            break
        fi
        total_time_us=$((total_time_us + time_output))
        all_times_us+="${time_output} "
        if [ "$min_time_us" -eq -1 ] || [ "$time_output" -lt "$min_time_us" ]; then min_time_us=$time_output; fi
        if [ "$time_output" -gt "$max_time_us" ]; then max_time_us=$time_output; fi
    done

    if [ "$run_status" == "PASS" ]; then
        avg_time_us=$((total_time_us / NUM_REPETITIONS))
        # Use awk for floating point division for milliseconds
        avg_time_ms=$(awk -v avg_us="$avg_time_us" 'BEGIN {printf "%.3f", avg_us / 1000}')
        min_time_ms=$(awk -v min_us="$min_time_us" 'BEGIN {printf "%.3f", min_us / 1000}')
        max_time_ms=$(awk -v max_us="$max_time_us" 'BEGIN {printf "%.3f", max_us / 1000}')
        printf "Result: PASS, Avg Time: %s ms (%d us) [Min: %s ms, Max: %s ms, Reps: %d]
Details (us): %s
" "$avg_time_ms" "$avg_time_us" "$min_time_ms" "$max_time_ms" "$NUM_REPETITIONS" "$all_times_us"
    else
        printf "Result: %s
" "$run_status"
        OVERALL_SERIAL_STATUS="FAIL"
    fi
done
echo "----------------------------"
echo "Overall Serial Test Status: ${OVERALL_SERIAL_STATUS}"
echo ""

echo "--- Running Parallel Tests (averaged over ${NUM_REPETITIONS} repetitions) ---"
OVERALL_PARALLEL_STATUS="PASS"
for threads in "${THREAD_COUNTS[@]}"; do
    echo "*** Using ${threads} OpenMP thread(s) ***"
    export OMP_NUM_THREADS=$threads
    STATUS_FOR_CURRENT_THREADS="PASS"

    for i in "${TEST_CASES[@]}"; do
        INPUT_FILE="${NTT_DATA_DIR}/${i}.in"
        OUTPUT_FILE="${NTT_DATA_DIR}/${i}.out"
        if [ ! -f "$INPUT_FILE" ] || [ ! -f "$OUTPUT_FILE" ]; then
            echo "Parallel Test Case ${i} (${INPUT_FILE##*/}, Threads: ${threads}): Warning: Files not found. Skipping."
            STATUS_FOR_CURRENT_THREADS="FAIL" # Or SKIP
            OVERALL_PARALLEL_STATUS="FAIL"
            continue
        fi

        echo -n "Parallel Test Case ${i} (${INPUT_FILE##*/}, Threads: ${threads}): "
        total_time_us=0
        run_status="PASS"
        min_time_us=-1
        max_time_us=0
        all_times_us=""

        for rep in $(seq 1 ${NUM_REPETITIONS}); do
            time_output=$("${PARALLEL_EXE}" "${INPUT_FILE}" "${OUTPUT_FILE}" 2>>"${ERROR_LOG_FILE}")
            exit_code=$?
            if [ $exit_code -ne 0 ]; then
                run_status="FAIL (run ${rep}, code ${exit_code})"
                echo "" # Newline
                echo "Error during Parallel Test Case ${i}, repetition ${rep}, Threads ${threads}. Check ${ERROR_LOG_FILE}" >> /dev/stderr
                break
            fi
            total_time_us=$((total_time_us + time_output))
            all_times_us+="${time_output} "
            if [ "$min_time_us" -eq -1 ] || [ "$time_output" -lt "$min_time_us" ]; then min_time_us=$time_output; fi
            if [ "$time_output" -gt "$max_time_us" ]; then max_time_us=$time_output; fi
        done

        if [ "$run_status" == "PASS" ]; then
            avg_time_us=$((total_time_us / NUM_REPETITIONS))
            avg_time_ms=$(awk -v avg_us="$avg_time_us" 'BEGIN {printf "%.3f", avg_us / 1000}')
            min_time_ms=$(awk -v min_us="$min_time_us" 'BEGIN {printf "%.3f", min_us / 1000}')
            max_time_ms=$(awk -v max_us="$max_time_us" 'BEGIN {printf "%.3f", max_us / 1000}')
            printf "Result: PASS, Avg Time: %s ms (%d us) [Min: %s ms, Max: %s ms, Reps: %d]
Details (us): %s
" "$avg_time_ms" "$avg_time_us" "$min_time_ms" "$max_time_ms" "$NUM_REPETITIONS" "$all_times_us"
        else
            printf "Result: %s
" "$run_status"
            STATUS_FOR_CURRENT_THREADS="FAIL"
            OVERALL_PARALLEL_STATUS="FAIL"
        fi
    done
    echo "Status for ${threads} threads: ${STATUS_FOR_CURRENT_THREADS}"
    echo "--------------------------------"
done

echo "--- Overall Test Summary ---"
if [ "$OVERALL_SERIAL_STATUS" == "PASS" ] && [ "$OVERALL_PARALLEL_STATUS" == "PASS" ]; then
    echo "Overall Status: ALL TESTS PASSED"
    # rm -f "${ERROR_LOG_FILE}" # Optionally remove error log on full success
    exit 0
else
    echo "Overall Status: SOME TESTS FAILED. Check ${ERROR_LOG_FILE} for details from C++ program stderr."
    exit 1
fi 