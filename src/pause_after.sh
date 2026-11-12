#!/usr/bin/env bash
# Pause the scaling run at the next clean rule boundary.
#
# Usage: pause_after.sh <target-file> <unit-name> [<log>]
#
# Waits until <target-file> exists, then stops the systemd user unit
# <unit-name>. Snakemake writes a rule's outputs only when the rule succeeds,
# so stopping right after the target appears never discards completed work.
set -u

TARGET="${1:?target file required}"
UNIT="${2:?unit name required}"
LOG="${3:-logs/pause.log}"

mkdir -p "$(dirname "$LOG")"
say() { echo "$(date -Is) $*" >> "$LOG"; }

say "watching for $TARGET (unit=$UNIT)"
while systemctl --user is-active --quiet "$UNIT"; do
    if [[ -e "$TARGET" ]]; then
        # Let Snakemake finish bookkeeping for the completed rule.
        sleep 8
        if systemctl --user is-active --quiet "$UNIT"; then
            systemctl --user stop "$UNIT"
            say "paused: $TARGET appeared, stopped $UNIT"
        else
            say "$UNIT already stopped after $TARGET appeared"
        fi
        exit 0
    fi
    sleep 10
done
say "$UNIT stopped before $TARGET appeared"
