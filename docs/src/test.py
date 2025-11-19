import re
with open("pyproject.toml") as f:
    text = f.read()
match = re.search(r'version\s*=\s*"([^"]+)"', text)
if match:
    print(match.group(1))
else:
    print("Version not found")
