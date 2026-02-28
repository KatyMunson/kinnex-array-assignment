#!/usr/bin/env bash
# =============================================================================
# commit_configs.sh — Config file git tracking helper
#
# Config files contain local paths and are normally hidden from git using
# --assume-unchanged. Use this script to temporarily expose them for committing
# template changes, or to restore the repo defaults.
#
# Usage:
#   bash commit_configs.sh --protect    # hide configs from git (default state)
#   bash commit_configs.sh --upload     # expose configs so they can be committed
#   bash commit_configs.sh --restore    # pull config defaults down from git
#
# Author:  KM
# Created: 2026-02
# =============================================================================

set -euo pipefail

# All config files that contain local paths and should normally be hidden
CONFIGS=(
    "prod_script/config.yaml"
    "test_script/config_qc.yaml"
    "test_script/config_subsample.yaml"
)

# Make sure we're running from the repo root
if [ ! -f ".gitignore" ]; then
    echo "ERROR: Run this script from the repo root (kinnex_array_assignment/)"
    exit 1
fi

usage() {
    echo "Usage: bash commit_configs.sh [--protect | --upload | --restore]"
    echo ""
    echo "  --protect   Hide config files from git (normal working state)"
    echo "  --upload    Expose config files so changes can be staged and committed"
    echo "  --restore   Discard local changes and pull config defaults from git"
    exit 1
}

if [ $# -ne 1 ]; then
    usage
fi

case "$1" in

    --protect)
        echo "Protecting config files from git tracking..."
        for f in "${CONFIGS[@]}"; do
            if [ -f "$f" ]; then
                git update-index --assume-unchanged "$f"
                echo "  protected: $f"
            else
                echo "  skipped (not found): $f"
            fi
        done
        echo "Done. Local changes to config files will not appear in git status."
        ;;

    --upload)
        echo "Exposing config files for git staging..."
        for f in "${CONFIGS[@]}"; do
            if [ -f "$f" ]; then
                git update-index --no-assume-unchanged "$f"
                echo "  unprotected: $f"
            else
                echo "  skipped (not found): $f"
            fi
        done
        echo ""
        echo "Config files are now visible to git."
        echo "Before committing, make sure they contain placeholder paths, not real ones."
        echo "Run 'git diff' to review changes, then 'git add' and 'git commit' as normal."
        echo "Run 'bash commit_configs.sh --protect' again when done."
        ;;

    --restore)
        echo "Restoring config files to repo defaults..."
        echo "WARNING: This will discard any local changes to config files."
        read -r -p "Are you sure? (y/N) " confirm
        if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
            echo "Aborted."
            exit 0
        fi
        for f in "${CONFIGS[@]}"; do
            # Must un-assume-unchanged before checkout will work
            git update-index --no-assume-unchanged "$f" 2>/dev/null || true
            if git ls-files --error-unmatch "$f" &>/dev/null; then
                git checkout HEAD -- "$f"
                # Re-protect after restoring
                git update-index --assume-unchanged "$f"
                echo "  restored: $f"
            else
                echo "  skipped (not in repo): $f"
            fi
        done
        echo "Done. Config files restored to repo defaults and re-protected."
        ;;

    *)
        usage
        ;;

esac
