# TRONKO BUILD/TEST COMMANDS

## Build Commands
- Build tronko-build: `cd tronko-build && make`
- Build tronko-assign: `cd tronko-assign && make`
- Debug builds: `make debug`
- Clean builds: `make clean`

## Examples
- Test tronko-build (single tree): `tronko-build -l -t example.tree -m example.fasta -x taxonomy.txt -d output_dir`
- Test tronko-assign: `tronko-assign -r -f database.txt -a reference.fasta -s -g reads.fasta -o results.txt -w`

# CODE STYLE GUIDELINES

## General
- Language: C (GNU99 standard)
- Indentation: 4 spaces (mixed with tabs in places)
- Memory management: Manual with malloc/free

## Naming Conventions
- Files: snake_case (e.g., readreference.c)
- Functions: camelCase (e.g., readReferenceTree)
- Constants/Macros: UPPERCASE (e.g., MAXQUERYLENGTH)
- Structs: lowercase (e.g., node, queryMatPaired)

## Error Handling
- Use fprintf(stderr) for error messages
- Exit with negative codes for critical failures
- Check return values for NULL pointers