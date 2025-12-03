#!/bin/bash
echo -e "electronID Cut for MC"
echo -n "Enter the Date: "
read date
echo -n "Attempt: "
read Attempt

# === Input directories ===
INPUT_DIR="/group/belle/users/amubarak/PID_FullRange/Ds_${date}25_${Attempt}_ccbar_ReverseID"
SUBDIRS=("sub00" "sub01")

# === Output directory ===
OUTPUT_DIR="/group/belle/users/amubarak/02-Grid/Sample_Grid_ReverseID/Ds_${date}25_${Attempt}_ccbar_ReverseID_Chunks"
mkdir -p "$OUTPUT_DIR"

# === Collect all root files ===
all_files=()
for sub in "${SUBDIRS[@]}"; do
    while IFS= read -r -d '' file; do
        all_files+=("$file")
    done < <(find "$INPUT_DIR/$sub" -name "*.root" -print0)
done

total_files=${#all_files[@]}
chunks=20
files_per_chunk=$(( (total_files + chunks - 1) / chunks ))  # ceil division

echo "🔍 Found $total_files ROOT files"
echo "🔧 Merging into $chunks chunks ($files_per_chunk files each)"

# === Merge using hadd ===
for ((i = 0; i < chunks; i++)); do
    start=$((i * files_per_chunk))
    end=$((start + files_per_chunk))
    chunk_files=("${all_files[@]:start:files_per_chunk}")
    
    output_file="$OUTPUT_DIR/ccbar_chunk_$(printf "%02d" $i).root"
    
    if [ ${#chunk_files[@]} -eq 0 ]; then
        echo "⚠️  Chunk $i has no files, skipping."
        continue
    fi

    echo "📦 Merging chunk $i with ${#chunk_files[@]} files into: $output_file"
    hadd -f "$output_file" "${chunk_files[@]}" > /dev/null 2>&1

    if [ $? -ne 0 ]; then
        echo "❌ Failed to create chunk $i"
    else
        echo "✅ Chunk $i written"
    fi
done

echo "🏁 All chunks completed."