#!/usr/bin/env bash
set -e

# Fix ownership of output files
finish() {
    # Fix ownership of output files
    user_id=$(stat -c '%u:%g' /data)
    chown -R ${user_id} /data
}
trap finish EXIT

# gather time and memory statistics, save everything to log
echo "/usr/bin/time -f '\\\\nDEBUG_MAX_MEM:%M\\\\nDEBUG_RUNTIME:%E\\\\n' /opt/shasta/build/shasta-install/bin/shasta $@\n" > /data/shasta.log
eval "/usr/bin/time -f '\\nDEBUG_MAX_MEM:%M\\nDEBUG_RUNTIME:%E\\n' /opt/shasta/build/shasta-install/bin/shasta $@" 2>&1 | tee -a /data/shasta.log
