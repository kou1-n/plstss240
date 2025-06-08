# Agent Memory Guidelines

## File Search Best Practices

### 1. Primary Search Method
- Always use `file_search` as the first method when looking for specific files
- This provides the most reliable way to locate files when the filename is known
- Example: `file_search` for "README.md" or other specific files

### 2. Directory Listing Considerations
- When using `list_dir`, be aware that:
  - With many files, some might be missed due to display limits
  - Consider using multiple search methods for verification
  - Be extra careful when searching for important files like README.md

### 3. Fallback Procedures
- If a file cannot be found:
  - Try multiple search methods
  - Consider different possible file names or locations
  - Ask for user confirmation if necessary

### 4. Verification Steps
- Always verify file existence through multiple methods when possible
- Double-check important files like documentation
- Maintain awareness of potential display limitations in directory listings

These guidelines help ensure more reliable file discovery and prevent missing important files in the codebase. 