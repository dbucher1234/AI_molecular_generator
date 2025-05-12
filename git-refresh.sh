#!/usr/bin/env bash
#
# git-refresh.sh  —  one-command add → commit → pull-rebase → push
# Usage: ./git-refresh.sh "my commit message"
# If no argument is given, the commit message defaults to "Quick sync".

set -e  # abort on first error

msg=${1:-"Quick sync"}                    # use arg1 or default message
branch=$(git symbolic-ref --short HEAD)   # current branch name

echo "🔄  Syncing branch: $branch → origin/$branch"

# Stage all changes
git add -A

# Commit if there is anything to commit
if git diff --cached --quiet; then
    echo "🛈  Nothing to commit."
else
    git commit -m "$msg"
fi

# If no upstream is set, push and set it
if ! git rev-parse --abbrev-ref --symbolic-full-name "@{u}" &>/dev/null; then
    echo "⤴️  Setting upstream to origin/$branch"
    git push -u origin "$branch"
else
    # Otherwise pull-rebase then push
    git pull --rebase origin "$branch"
    git push origin "$branch"
fi

echo "✅  Done."

