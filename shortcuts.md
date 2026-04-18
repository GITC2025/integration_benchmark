# shortcuts

## md5check shortcut
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'

md5check() {
  local local_file="${1:-localchecksumfile.md5}"
  local target_file="${2:-targetchecksumfile.md5}"

  awk '
  NR==FNR {
    ena[$2] = $1
    next
  }
  {
    loc = $1
    $1 = ""
    sub(/^ +/, "")
    f = $0
    print f, (ena[f] == loc ? "MATCH" : "MISMATCH")
  }
  ' "$target_file" "$local_file"
}
EOF
source ~/.bashrc
```
```bash
# use shortcut
md5check localchecksumfile.md5 targetchecksumfile.md5
```

## unzipfastq and move fastq.txt to current folder shortcut
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'

unzipfastq() {
  # unzip all zip files in current directory
  for z in *.zip; do
    [ -f "$z" ] || continue
    unzip -o "$z"
  done

  # move each fastqc_data.txt from subdirectories into current directory
  # prefix with parent directory name to avoid collisions
  find . -mindepth 2 -type f -name 'fastqc_data.txt' -print0 \
    | while IFS= read -r -d '' f; do
        d=$(basename "$(dirname "$f")")
        mv -n "$f" "./${d}_fastqc_data.txt"
      done
}
EOF
source ~/.bashrc
```


## FastQC txt parser shortcut
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'
fastqcfail() {
  shopt -s nullglob
  local files=( *.txt )
  [ ${#files[@]} -gt 0 ] || { echo "no *.txt files here"; return 1; }

  awk 'BEGIN{OFS="\t"; print "FILE","MODULE","STATUS"}
  /^FAIL[[:space:]]/{
    status=$1; $1=""; sub(/^[[:space:]]+/,""); print FILENAME, $0, status
  }
  /^>>/{
    line=$0; sub(/^>>/,"",line); split(line,a,"\t");
    if(toupper(a[2])=="FAIL") print FILENAME, a[1], a[2]
  }' "${files[@]}"
}
EOF
source ~/.bashrc
```
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'
fastqcwarn() {
  shopt -s nullglob
  local files=( *.txt )
  [ ${#files[@]} -gt 0 ] || { echo "no *.txt files here"; return 1; }

  awk 'BEGIN{OFS="\t"; print "FILE","MODULE","STATUS"}
  /^WARN[[:space:]]/{
    status=$1; $1=""; sub(/^[[:space:]]+/,""); print FILENAME, $0, status
  }
  /^>>/{
    line=$0; sub(/^>>/,"",line); split(line,a,"\t");
    if(toupper(a[2])=="WARN") print FILENAME, a[1], a[2]
  }' "${files[@]}"
}
EOF
source ~/.bashrc
```
```bash
# use shortcuts
fastqcfail
fastqcwarn
```
## zcat check chemistry
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'
zcatR1() {
  local files=("$@")
  if [ ${#files[@]} -eq 0 ]; then files=( *1.fastq.gz ); fi
  for f in "${files[@]}"; do
    echo "$f:"
    zcat "$f" | head -n 10 | awk 'NR%4==2 {print length($0)}' | sort -u
  done
}
EOF
source ~/.bashrc
```
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'
zcatR2() {
  local files=("$@")
  if [ ${#files[@]} -eq 0 ]; then files=( *2.fastq.gz ); fi
  for f in "${files[@]}"; do
    echo "$f:"
    zcat "$f" | head -n 10 | awk 'NR%4==2 {print length($0)}' | sort -u
  done
}
EOF
source ~/.bashrc
```
## cellranger QC parser for scRNA data
```bash
under construction
```

## cellranger pipestance completion parser
```bash
printf '\n' >> ~/.bashrc
cat << 'EOF' >> ~/.bashrc

pipesuccess() {
local failed_files=()
local search_path="${1:-.}"
[[ "$search_path" != */ ]] && search_path="$search_path/"
for log in "${search_path}"*.out; do
if ! tail -n 5 "$log" | grep -q "Pipestance completed successfully!"; then
failed_files+=("$log")
fi
done
if [ ${#failed_files[@]} -eq 0 ]; then
echo "All pipestance completed successfully."
else
echo "Warning: The following jobs did NOT complete successfully:"
printf '%s\n' "${failed_files[@]}"
return 1
fi
}
EOF
source ~/.bashrc
```

## other shortcuts
```bash
alias myjobs='queue -u "$USER'
alias smallnode='salloc --time=04:00:00 --cpus-per-task=8 --mem=16G'
alias bignode='salloc --time=04:00:00 --cpus-per-task=8 --mem=128G'
# dataset 1
alias GSE='cd /global/scratch/$USER/12marchENA_GSE'

# dataset 2
alias SRP='cd /global/scratch/$USER/4samplesSRP351272'

# bind the alias to persist across sessions
echo "alias myalias='command here'" >> ~/.bashrc
source ~/.bashrc
```

## bind shell fxns
```bash
printf '\n' >> ~/.bashrc
cat >> ~/.bashrc <<'EOF'

myfunc() {
  # function body
  echo "doing something"
}
EOF
source ~/.bashrc
```




