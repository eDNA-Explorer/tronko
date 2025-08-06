#!/bin/bash

echo "=== Testing tronko-assign with Example Data and Logging ==="
echo ""

cd /app

# First, let's check if the example data exists
echo "1. Checking for example data files..."
if [ ! -f "tronko-build/example_datasets/single_tree/reference_tree.txt" ]; then
    echo "❌ Reference tree file not found. Need to build the database first."
    echo "Building the single tree example database..."
    cd tronko-build
    ./tronko-build -l \
        -m example_datasets/single_tree/Charadriiformes_MSA.fasta \
        -x example_datasets/single_tree/Charadriiformes_taxonomy.txt \
        -t example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
        -d example_datasets/single_tree
    cd ..
fi

if [ ! -f "example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta" ]; then
    echo "❌ Test reads file not found!"
    exit 1
fi

echo "✅ Example data files found!"
echo ""

# Test 1: Basic run without logging (baseline)
echo "2. Testing basic tronko-assign without logging (baseline)..."
cd tronko-assign
./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s \
    -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o test_results_baseline.txt \
    -w

if [ $? -eq 0 ]; then
    echo "✅ Baseline test successful!"
    echo "Results written to test_results_baseline.txt"
    wc -l test_results_baseline.txt
else
    echo "❌ Baseline test failed!"
    exit 1
fi
echo ""

# Test 2: Run with INFO level logging
echo "3. Testing with INFO level logging (-V 2)..."
./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s \
    -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o test_results_info.txt \
    -w \
    -V 2

if [ $? -eq 0 ]; then
    echo "✅ INFO logging test successful!"
    echo "Results written to test_results_info.txt"
    wc -l test_results_info.txt
else
    echo "❌ INFO logging test failed!"
fi
echo ""

# Test 3: Run with DEBUG level logging and timing
echo "4. Testing with DEBUG level logging, timing, and log file (-V 3 -T -l)..."
./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s \
    -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o test_results_debug.txt \
    -w \
    -V 3 \
    -T \
    -l debug_timing.log

if [ $? -eq 0 ]; then
    echo "✅ DEBUG logging with timing test successful!"
    echo "Results written to test_results_debug.txt"
    echo "Log written to debug_timing.log"
    echo ""
    echo "=== Debug Log Contents ==="
    cat debug_timing.log
    echo "=========================="
    wc -l test_results_debug.txt
else
    echo "❌ DEBUG logging test failed!"
fi
echo ""

# Test 4: Run with resource monitoring
echo "5. Testing with resource monitoring (-R)..."
./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s \
    -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o test_results_resources.txt \
    -w \
    -V 2 \
    -R \
    -l resource_monitor.log

if [ $? -eq 0 ]; then
    echo "✅ Resource monitoring test successful!"
    echo "Results written to test_results_resources.txt"
    echo "Resource log written to resource_monitor.log"
    echo ""
    echo "=== Resource Monitor Log Contents ==="
    cat resource_monitor.log
    echo "====================================="
    wc -l test_results_resources.txt
else
    echo "❌ Resource monitoring test failed!"
fi
echo ""

# Test 5: Full logging features test
echo "5. Testing all logging features together (-V 3 -R -T -l)..."
./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s \
    -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o test_results_full.txt \
    -w \
    -V 3 \
    -R \
    -T \
    -l full_logging.log

if [ $? -eq 0 ]; then
    echo "✅ Full logging features test successful!"
    echo "Results written to test_results_full.txt"
    echo "Full log written to full_logging.log"
    echo ""
    echo "=== Full Logging Log Contents (last 50 lines) ==="
    tail -50 full_logging.log
    echo "=================================================="
    wc -l test_results_full.txt
else
    echo "❌ Full logging test failed!"
fi
echo ""

# Test 6: Crash testing with corrupted paired-end data
echo "6. Testing crash debugging with corrupted paired-end data..."
echo "This test is expected to potentially crash to verify our crash debugging system works!"

# Check if the paired-end files exist
if [ ! -f "../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta" ]; then
    echo "❌ Corrupted paired-end read1 file not found!"
else
    echo "✅ Found corrupted paired-end read1 file"
    
    # Test with corrupted data and full crash debugging enabled
    echo "Running with corrupted data and full crash debugging..."
    ./tronko-assign -r \
        -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
        -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
        -p \
        -1 ../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta \
        -2 ../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta \
        -o test_results_crash_test.txt \
        -w \
        -V 3 \
        -R \
        -T \
        -l crash_test.log
    
    crash_exit_code=$?
    if [ $crash_exit_code -eq 0 ]; then
        echo "✅ Corrupted data test completed successfully (no crash occurred)"
        echo "Results written to test_results_crash_test.txt"
        wc -l test_results_crash_test.txt
    else
        echo "💥 Process crashed with exit code $crash_exit_code (this was expected!)"
        echo "Checking for crash reports..."
        ls -la /tmp/tronko_assign_crash_*.crash 2>/dev/null || echo "No crash reports found"
        
        if [ -f "crash_test.log" ]; then
            echo ""
            echo "=== Crash Test Log Contents (last 20 lines) ==="
            tail -20 crash_test.log
            echo "==============================================="
        fi
    fi
fi
echo ""

# Verify all results are consistent
echo "8. Verifying result consistency across different logging levels..."
if diff test_results_baseline.txt test_results_info.txt > /dev/null; then
    echo "✅ Baseline and INFO results are identical"
else
    echo "❌ Baseline and INFO results differ!"
fi

if diff test_results_baseline.txt test_results_debug.txt > /dev/null; then
    echo "✅ Baseline and DEBUG results are identical"
else
    echo "❌ Baseline and DEBUG results differ!"
fi

if diff test_results_baseline.txt test_results_full.txt > /dev/null; then
    echo "✅ All results are consistent - logging doesn't affect output!"
else
    echo "❌ Results differ between logging modes!"
fi

echo ""
echo "=== Test Summary ==="
echo "All tests completed. Check the log files for detailed logging output:"
ls -la *.log *.txt
echo ""
echo "=== Logging Implementation Test Complete ==="