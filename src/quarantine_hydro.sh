#!/bin/bash

# Hydro file quarantine script for 6-hour timestep testing
# Moves files to quarantine directory (does NOT delete)

# Directories
SOURCE_DIR="/ddnA/project/holstein/NESTS/CMS_test_nests"
QUARANTINE_DIR="/ddnA/project/holstein/NESTS/CMS_test_nests/quarantine"

# Create quarantine directory if it doesn't exist
mkdir -p "$QUARANTINE_DIR"

# Set mode: "test" for debugging, "full" for March 1-17 2019
MODE="test"  # Change to "full" when ready

echo "=== Hydro File Quarantine Script ==="
echo "Mode: $MODE"
echo "Source: $SOURCE_DIR"
echo "Quarantine: $QUARANTINE_DIR"
echo ""

# Hours to KEEP for 6-hour timesteps (all others will be moved to quarantine)
KEEP_HOURS=("00" "06" "12" "18")

# Function to check if a given hour should be kept
should_keep_hour() {
    local hour="$1"
    for keep_hour in "${KEEP_HOURS[@]}"; do
        if [[ "$hour" == "$keep_hour" ]]; then
            return 0  # Keep this hour
        fi
    done
    return 1  # Quarantine this hour
}

# Function to extract hour from filename
# For filename like "nest_1_20190301140000.nc", hour "14" is at positions 15-16
extract_hour() {
    local filename="$1"
    echo "${filename:15:2}"
}

if [[ "$MODE" == "test" ]]; then
    echo "=== TEST MODE: Processing sample files for debugging ==="
    
    # Test with a few specific files to verify logic
    TEST_PATTERNS=(
        "nest_1_20190301000*.nc"  # Should KEEP (hour 00)
        "nest_1_20190301010*.nc"  # Should MOVE (hour 01)
        "nest_1_20190301060*.nc"  # Should KEEP (hour 06)
        "nest_1_20190301070*.nc"  # Should MOVE (hour 07)
        "nest_1_20190301120*.nc"  # Should KEEP (hour 12)
    )
    
    echo "Testing hour extraction and logic..."
    files_found=0
    files_moved=0
    files_kept=0
    
    for pattern in "${TEST_PATTERNS[@]}"; do
        echo "Checking pattern: $pattern"
        
        for file in $SOURCE_DIR/$pattern; do
            if [[ -f "$file" ]]; then
                files_found=$((files_found + 1))
                filename=$(basename "$file")
                
                # Extract hour from filename
                hour=$(extract_hour "$filename")
                
                echo "  File: $filename"
                echo "    Extracted hour: '$hour'"
                
                if should_keep_hour "$hour"; then
                    files_kept=$((files_kept + 1))
                    echo "    Decision: KEEP (6-hour timestep)"
                else
                    echo "    Decision: MOVE to quarantine"
                    mv "$file" "$QUARANTINE_DIR/"
                    if [[ $? -eq 0 ]]; then
                        files_moved=$((files_moved + 1))
                        echo "    Status: ✓ Moved successfully"
                    else
                        echo "    Status: ✗ Error moving file"
                    fi
                fi
                echo ""
            fi
        done
    done
    
    echo "TEST SUMMARY:"
    echo "Files found: $files_found"
    echo "Files kept: $files_kept" 
    echo "Files moved to quarantine: $files_moved"
    echo "Hours that should be kept: ${KEEP_HOURS[*]}"
    
elif [[ "$MODE" == "full" ]]; then
    echo "=== FULL MODE: Processing March 1-17, 2019 ==="
    
    # Process all days from March 1-17, 2019
    total_files_found=0
    total_files_moved=0
    total_files_kept=0
    
    for day in {01..17}; do
        echo "Processing March $day, 2019..."
        
        pattern="nest_1_201903${day}*.nc"
        day_files_found=0
        day_files_moved=0
        day_files_kept=0
        
        for file in $SOURCE_DIR/$pattern; do
            if [[ -f "$file" ]]; then
                total_files_found=$((total_files_found + 1))
                day_files_found=$((day_files_found + 1))
                filename=$(basename "$file")
                
                # Extract hour from filename
                hour=$(extract_hour "$filename")
                
                if should_keep_hour "$hour"; then
                    total_files_kept=$((total_files_kept + 1))
                    day_files_kept=$((day_files_kept + 1))
                    echo "  KEEP: $filename (hour: $hour)"
                else
                    echo "  MOVE: $filename (hour: $hour)"
                    mv "$file" "$QUARANTINE_DIR/"
                    if [[ $? -eq 0 ]]; then
                        total_files_moved=$((total_files_moved + 1))
                        day_files_moved=$((day_files_moved + 1))
                    else
                        echo "    ✗ Error moving $filename"
                    fi
                fi
            fi
        done
        
        echo "  Day $day summary: $day_files_found total, $day_files_kept kept, $day_files_moved moved"
        echo ""
    done
    
    echo "FULL PROCESSING SUMMARY:"
    echo "Total files processed: $total_files_found"
    echo "Files kept for 6-hour timesteps: $total_files_kept"
    echo "Files moved to quarantine: $total_files_moved"
    echo "Hours kept: ${KEEP_HOURS[*]}"
    
else
    echo "ERROR: Invalid mode '$MODE'"
    echo "Set MODE to 'test' or 'full' at the top of this script"
    exit 1
fi

echo ""
echo "=== Operation Complete ==="
echo "Quarantine directory: $QUARANTINE_DIR"

# Show a few examples of what's in quarantine
if [[ -d "$QUARANTINE_DIR" ]]; then
    quarantine_count=$(ls -1 "$QUARANTINE_DIR"/*.nc 2>/dev/null | wc -l)
    if [[ $quarantine_count -gt 0 ]]; then
        echo "Files in quarantine: $quarantine_count"
        echo "Sample quarantined files:"
        ls "$QUARANTINE_DIR"/*.nc | head -3
    fi
fi