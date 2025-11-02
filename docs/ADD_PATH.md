# Adding a Directory to PATH (Only if Not Already Present)

## Using pattern matching (cleanest)
```bash
if [[ ":$PATH:" != *":/new/directory:"* ]]; then
    export PATH="/new/directory:$PATH"
fi
```

## Using case statement
```bash
case ":$PATH:" in
    *:/new/directory:*)
        ;;
    *)
        export PATH="/new/directory:$PATH"
        ;;
esac
```

## As a reusable function (recommended for `.bashrc`)
```bash
add_to_path() {
    if [[ ":$PATH:" != *":$1:"* ]]; then
        export PATH="$1:$PATH"
    fi
}

# Usage:
add_to_path "/new/directory"
add_to_path "$HOME/bin"
```

## One-liner (compact)
```bash
[[ ":$PATH:" != *":/new/directory:"* ]] && export PATH="/new/directory:$PATH"
```

## Key points
- The colons (`:`) around `$PATH` and the directory prevent partial matches (e.g., `/usr/local/bin` matching `/usr/local`)
- Add to the beginning with `"$dir:$PATH"` or end with `"$PATH:$dir"`
- Quote variables to handle spaces in paths

The function approach is best for `.bashrc` since you can reuse it for multiple directories.
