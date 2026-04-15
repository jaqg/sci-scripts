#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -p PARAMETER [-v VALUE] [-d DIRECTORY] [-x EXCLUDE] [-c] [--delete] [-h]"
    echo "  -p PARAMETER    Name of the parameter to modify or delete (e.g., PROPERTY, MODEL)."
    echo "  -v VALUE        New value to assign to the parameter (ignored if --delete is used)."
    echo "  -d DIRECTORY    Directory pattern to search (default: current directory)."
    echo "  -x EXCLUDE      Comma-separated list of directory patterns to exclude."
    echo "  -c              Create fcc.inp if it does not exist (using 'fcclasses3 -h')."
    echo "  --delete        Delete the parameter line from the input file."
    echo "  -h              Display this help message."
    exit 1
}

# Default values
directory="."
excluded_dirs=()
create_if_missing=false
delete_param=false

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -p)
            parameter="$2"
            shift 2
            ;;
        -v)
            value="$2"
            shift 2
            ;;
        -d)
            directory="$2"
            shift 2
            ;;
        -x)
            IFS=',' read -r -a excluded_dirs <<< "$2"
            shift 2
            ;;
        -c)
            create_if_missing=true
            shift
            ;;
        --delete)
            delete_param=true
            shift
            ;;
        -h)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$parameter" ]]; then
    echo "Error: Missing required parameter argument."
    usage
fi

if ! $delete_param && [[ -z "$value" ]]; then
    echo "Error: -v VALUE is required when not using --delete."
    usage
fi

# Get directories matching the pattern using glob (safer than parsing ls output)
shopt -s nullglob
directories=($directory)
shopt -u nullglob

# Exit if no directories are found
if [[ ${#directories[@]} -eq 0 ]]; then
    echo "Error: No directories found matching pattern '$directory'"
    exit 1
fi

echo "Processing directories: ${directories[*]}"

# Escape special characters in value for sed (delimiter is #, so escape & and #).
# Note: values containing newlines are not supported and would break the substitution.
escaped_value=$(printf "%s" "$value" | sed 's|[&#]|\\&|g')

# Escape parameter for use as a BRE pattern in grep/sed.
escaped_param=$(printf '%s' "$parameter" | sed 's/[][\\.*^$]/\\&/g')

# Build prune args for excluded directories
prune_args=()
for dir in "${excluded_dirs[@]}"; do
    shopt -s nullglob
    expanded=($dir)
    shopt -u nullglob
    for expanded_dir in "${expanded[@]}"; do
        prune_args+=(-path "$expanded_dir" -prune -o)
    done
done

# Build the final type clause separately from prune args
if $create_if_missing; then
    type_args=(-type d -print)
else
    type_args=(-type f -name "fcc.inp" -print)
fi

# Execute `find` command to get the list of directories/files
find "${directories[@]}" "${prune_args[@]}" "${type_args[@]}" | while read -r dir; do
    if $create_if_missing; then
        file="$dir/fcc.inp"
        if [[ ! -f "$file" ]]; then
            echo "Creating fcc.inp in $dir..."
            (cd "$dir" && fcclasses3 -h > fcc.inp)
        fi
    else
        file="$dir"
    fi

    if $delete_param; then
        # Delete the parameter line if --delete is set
        if grep -q "^$escaped_param" "$file"; then
            sed -i "/^$escaped_param[[:space:]]*=.*/d" "$file"
            echo "Deleted $parameter from $file"
        else
            echo "$parameter not found in $file"
        fi
    else
        # Modify or add the parameter
        if grep -q "^$escaped_param" "$file"; then
            sed -i "s#^$escaped_param[[:space:]]*=.*#$parameter =   $escaped_value#" "$file"
            echo "Updated $parameter in $file"
        else
            echo "$parameter =   $value" >> "$file"
            echo "Added $parameter to $file"
        fi
    fi

done

echo "All modifications completed."
