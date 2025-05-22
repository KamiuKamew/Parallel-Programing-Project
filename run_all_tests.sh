#!/bin/bash

MAIN_CC_PATH="ntt/main.cc"
MAIN_CC_BAK_PATH="ntt/main.cc.bak"
OUTPUT_DIR="ntt/test_outputs"

# Ensure ntt directory exists
if [ ! -d "ntt" ]; then
  echo "Error: ntt directory not found at $(pwd)/ntt. Please run this script from the workspace root."
  exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"
echo "Output files will be saved in $(pwd)/${OUTPUT_DIR}"

# Backup main.cc
echo "Backing up ${MAIN_CC_PATH} to ${MAIN_CC_BAK_PATH}..."
if [ -f "${MAIN_CC_BAK_PATH}" ]; then
    echo "Warning: Backup file ${MAIN_CC_BAK_PATH} already exists. It will be overwritten."
fi
cp "${MAIN_CC_PATH}" "${MAIN_CC_BAK_PATH}"

# Trap to restore main.cc on exit, interrupt, or termination
cleanup() {
    echo ""
    echo "Cleaning up..."
    if [ -f "${MAIN_CC_BAK_PATH}" ]; then
        echo "Restoring ${MAIN_CC_PATH} from ${MAIN_CC_BAK_PATH}..."
        cp -f "${MAIN_CC_BAK_PATH}" "${MAIN_CC_PATH}"
        rm -f "${MAIN_CC_BAK_PATH}"
        echo "Backup file ${MAIN_CC_BAK_PATH} removed."
    else
        echo "Warning: Backup file ${MAIN_CC_BAK_PATH} not found. Cannot restore ${MAIN_CC_PATH}."
    fi
    echo "Cleanup complete."
}
trap cleanup EXIT INT TERM

cores_threads_configs=("1 1" "2 2" "4 4" "8 8")

u32_functions=("poly_multiply_ntt" "poly_multiply_ntt_omp" "poly_multiply_ntt_pthread_simple")
u64_functions=("poly_multiply_ntt_simd" "poly_multiply_ntt_omp" "poly_multiply_ntt_crt" "poly_multiply_ntt_pthread_crt" "poly_multiply_ntt_pthread_simple")

# --- Helper function to modify main.cc ---
# Arg1: datatype (e.g., u32)
# Arg2: function_to_enable (e.g., poly_multiply_ntt)
modify_main_cc() {
  local dtype="$1"
  local func_to_enable="$2"

  echo "  Modifying ${MAIN_CC_PATH} for ${dtype} and ${func_to_enable}"

  # Restore from original backup to ensure a clean state for sed commands
  cp "${MAIN_CC_BAK_PATH}" "${MAIN_CC_PATH}"

  # 1. Set _main<datatype> in the main() function (around line 147)
  sed -i "147s/_main<u[0-9a-zA-Z_]*>(argc, argv);/_main<${dtype}>(argc, argv);/" "${MAIN_CC_PATH}"
  if [ $? -ne 0 ]; then echo "Error setting datatype ${dtype}."; exit 1; fi

  # 2. Comment out the default active line (poly_multiply_ntt on line 113) from the backup's state.
  # This makes all poly_multiply_ lines (113-118) commented initially for the next step.
  # Assumes main.cc.bak has poly_multiply_ntt on line 113 uncommented and others on 114-118 commented.
  sed -i '113s|^    poly_multiply_ntt(a, b, ab, n_, p_);|    // poly_multiply_ntt(a, b, ab, n_, p_);|' "${MAIN_CC_PATH}"
  if [ $? -ne 0 ]; then echo "Error commenting default poly_multiply_ntt."; exit 1; fi
  
  # 3. Uncomment the target function (within lines 113-118)
  # The pattern searches for "    // FUNCTION_NAME(a, b, ab, n_, p_);" and removes the "// "
  # It must match the exact function name.
  local sed_command="113,118s|^    // \(${func_to_enable}(a, b, ab, n_, p_)\);|    \1;|"
  sed -i "${sed_command}" "${MAIN_CC_PATH}"
  if [ $? -ne 0 ]; then echo "Error enabling function ${func_to_enable}."; exit 1; fi

  # Verification (optional, for debugging script development)
  # echo "    Verification for ${dtype}, ${func_to_enable}:"
  # echo "    Active _main call:"
  # grep -n "_main<.*>(argc, argv);" "${MAIN_CC_PATH}"
  # echo "    Active poly_multiply function (lines 113-118):"
  # sed -n '113,118p' "${MAIN_CC_PATH}" | grep -v "//" | grep "poly_multiply_"
}
# --- End of helper function ---

echo "Starting test runs..."

for ct_config in "${cores_threads_configs[@]}"; do
  read -r cores threads <<< "${ct_config}"
  echo ""
  echo "===== Running tests for Cores: ${cores}, Threads: ${threads} ====="

  # Test u32 functions
  echo "  --- Processing u32 functions ---"
  for func_name in "${u32_functions[@]}"; do
    echo "    Testing u32, ${func_name}"
    modify_main_cc "u32" "${func_name}"

    echo "    Compiling..."
    (cd ntt && g++ main.cc -o main -O2 -fopenmp -lpthread -std=c++11)
    compile_status=$?
    if [ ${compile_status} -ne 0 ]; then
      echo "    Compilation FAILED for u32, ${func_name} with cores ${cores}, threads ${threads}. Skipping test."
      # Restore main.cc to a clean state for the next iteration just in case
      cp "${MAIN_CC_BAK_PATH}" "${MAIN_CC_PATH}"
      continue
    fi

    echo "    Running test.sh (param1=2, cores=${cores}, threads=${threads})..."
    (cd ntt && ./test.sh 2 "${cores}" "${threads}") # test.sh outputs to ntt/test.o
    test_status=$?
    if [ ${test_status} -ne 0 ]; then
      echo "    test.sh FAILED for u32, ${func_name} with cores ${cores}, threads ${threads}."
    else
      output_filename="test_cores${cores}_threads${threads}_u32_${func_name}.o"
      echo "    Saving output to ${OUTPUT_DIR}/${output_filename}"
      mv "ntt/test.o" "${OUTPUT_DIR}/${output_filename}"
    fi
    echo "    ------------------------------------"
  done

  # Test u64 functions
  echo "  --- Processing u64 functions ---"
  for func_name in "${u64_functions[@]}"; do
    echo "    Testing u64, ${func_name}"
    modify_main_cc "u64" "${func_name}"

    echo "    Compiling..."
    (cd ntt && g++ main.cc -o main -O2 -fopenmp -lpthread -std=c++11)
    compile_status=$?
    if [ ${compile_status} -ne 0 ]; then
      echo "    Compilation FAILED for u64, ${func_name} with cores ${cores}, threads ${threads}. Skipping test."
      cp "${MAIN_CC_BAK_PATH}" "${MAIN_CC_PATH}"
      continue
    fi

    echo "    Running test.sh (param1=2, cores=${cores}, threads=${threads})..."
    (cd ntt && ./test.sh 2 "${cores}" "${threads}")
    test_status=$?
    if [ ${test_status} -ne 0 ]; then
      echo "    test.sh FAILED for u64, ${func_name} with cores ${cores}, threads ${threads}."
    else
      output_filename="test_cores${cores}_threads${threads}_u64_${func_name}.o"
      echo "    Saving output to ${OUTPUT_DIR}/${output_filename}"
      mv "ntt/test.o" "${OUTPUT_DIR}/${output_filename}"
    fi
    echo "    ------------------------------------"
  done
done

echo ""
echo "All planned tests initiated."
echo "Output files should be in $(pwd)/${OUTPUT_DIR}"
echo "The cleanup trap will now restore ${MAIN_CC_PATH}." 