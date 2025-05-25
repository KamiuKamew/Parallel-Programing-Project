#!/bin/bash

cd ntt
g++ main.cc -o main -O2 -fopenmp -lpthread -std=c++11
if [ $? -ne 0 ]; then
  echo "编译失败，跳过执行 test.sh"
  exit 1
fi

RESULTS_FILE="test_results.txt"
# Clear previous results or create the file
echo "Test run started at $(date)" > "$RESULTS_FILE"
echo "Results will be saved to $RESULTS_FILE"
echo "" >> "$RESULTS_FILE" # Add a blank line for readability

# Loop for core/thread counts 1 to 8
for i in {1..8}
do
  # Output to console and file
  echo "Processing for test.sh 2 $i $i" | tee -a "$RESULTS_FILE"

  # Initialize arrays for sums and counts for the 5 test points
  declare -a sum_lat=(0 0 0 0 0)
  declare -a count_lat=(0 0 0 0 0)
  # Arrays to store N and P values for each test point (captured from the first successful run)
  declare -a test_point_n_values=()
  declare -a test_point_p_values=()

  # Repeat the test 5 times
  for j in {1..5}
  do
    echo "  Running repetition $((j))/5 for 2 $i $i..." # Keep this to console only for progress
    output=$(./test.sh 2 $i $i)
    
    num_awk_lines=0

    # Process awk output within a construct that keeps variable scope
    while IFS= read -r awk_line; do
        if [ -n "$awk_line" ]; then # Ensure line is not empty
            n_val=$(echo "$awk_line" | awk '{print $1}')
            p_val=$(echo "$awk_line" | awk '{print $2}')
            lat_val=$(echo "$awk_line" | awk '{print $3}')
            p_idx=$(echo "$awk_line" | awk '{print $4}') # 0-indexed point_idx

            if [ -z "$n_val" ] || [ -z "$p_val" ] || [ -z "$lat_val" ] || [ -z "$p_idx" ]; then
                # Output warning to console and file
                echo "    Warning: Malformed line from awk for repetition $j, params 2 $i $i: '$awk_line'" | tee -a "$RESULTS_FILE"
                continue 
            fi

            # Store N and P values if not already stored for this point_idx
            if [ -z "${test_point_n_values[$p_idx]}" ]; then
                test_point_n_values[$p_idx]=$n_val
                test_point_p_values[$p_idx]=$p_val
            fi
            
            current_sum=${sum_lat[$p_idx]}
            sum_lat[$p_idx]=$(echo "$current_sum + $lat_val" | bc -l)
            current_count=${count_lat[$p_idx]}
            count_lat[$p_idx]=$((current_count + 1))
            
            num_awk_lines=$((num_awk_lines + 1))
        fi
    done < <(echo "$output" | awk '''
        BEGIN { point_idx = 0; }
        /average latency for n = .* p = .* : .* \(us\)/ {
            if (point_idx < 5) {
                match($0, /n = ([0-9]+)[[:space:]]+p = ([0-9]+)[[:space:]]+: ([0-9\.]+)[[:space:]]+\(us\)/, arr)
                if (arr[1] && arr[2] && arr[3]) { 
                    printf "%s %s %s %s\n", arr[1], arr[2], arr[3], point_idx;
                }
                point_idx++; 
            }
        }
    ''')

    if [ "$num_awk_lines" -ne 5 ]; then
        # Output warning to console and file
        echo "    Warning: Repetition $j for 2 $i $i: Expected 5 test points, awk extracted $num_awk_lines." | tee -a "$RESULTS_FILE"
    fi
  done # End of j loop (repetitions)

  # Calculate and print averages
  echo "  Average latencies for test.sh 2 $i $i (over 5 repetitions):" | tee -a "$RESULTS_FILE"
  all_points_averaged_successfully=true
  for k in {0..4}
  do
    # Provide default "N/A" if N/P values were not captured for any reason
    n_k=${test_point_n_values[$k]:-N/A}
    p_k=${test_point_p_values[$k]:-N/A}
    
    if [ "${count_lat[$k]}" -gt 0 ]; then
      avg_lat=$(echo "scale=5; ${sum_lat[$k]} / ${count_lat[$k]}" | bc -l)
      # Output to console and file
      echo "    Test Point $((k+1)) (n=$n_k, p=$p_k): $avg_lat (us) (from ${count_lat[$k]} successful extractions)" | tee -a "$RESULTS_FILE"
    else
      # Output to console and file
      echo "    Test Point $((k+1)) (n=$n_k, p=$p_k): Error calculating average (0 valid extractions found)." | tee -a "$RESULTS_FILE"
      all_points_averaged_successfully=false
    fi
  done
  
  if ! $all_points_averaged_successfully; then
      # Output to console and file
      echo "    Note: Some test points may not have averages due to issues in data extraction from test.sh output." | tee -a "$RESULTS_FILE"
  fi
  echo "" | tee -a "$RESULTS_FILE" # Add a blank line for readability in the file
done # End of i loop (core/thread counts)

echo "All tests finished. Results saved to $RESULTS_FILE"
echo "All tests finished. Results appended to $RESULTS_FILE" >> "$RESULTS_FILE"
