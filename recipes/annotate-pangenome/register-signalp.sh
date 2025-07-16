#!/bin/bash

set -e

# === HELP MENU ===
function print_help() {
    echo ""
    echo "Usage: $(basename "$0") <archive.tar.gz>"
    echo ""
    echo "This script registers SignalP 4.1 from a provided archive and patches it for use in workflows."
    echo ""
    echo "Arguments:"
    echo "  <archive.tar.gz>   Path to the SignalP 4.1 archive file"
    echo ""
    echo "Example:"
    echo "  register-signalp.sh signalp-4.1g.Linux.tar.gz"
    echo ""
    echo "After installation, you can include SignalP in your Nextflow run with:"
    echo "  --signalp_path \"/your/path/to/signalp-4.1g-work/signalp-4.1g-install/\""
    echo ""
}

# === CHECK ARGUMENTS ===
if [[ "$1" == "-h" || "$1" == "--help" || -z "$1" ]]; then
    print_help
    exit 0
fi

ARCHIVE=$1
WORK_DIR=$(pwd)/"signalp-4.1g-work"
TARGET_DIR="signalp-4.1g-install"

# === USER CONFIGURATION ===
EXTRACTED_DIR_CALLED="$(basename $(tar -tf "${ARCHIVE}" | head -n 1))"

if [ -d "$WORK_DIR" ]; then
    echo "Directory $WORK_DIR exists. Deleting..."
    rm -rf "${WORK_DIR}"
    echo "Directory deleted."
else
    echo "Directory $WORK_DIR does not exist. Nothing to delete."
fi

mkdir -p "${WORK_DIR}"
tar --no-same-owner --directory=${WORK_DIR} -zxf "${ARCHIVE}"
cd "${WORK_DIR}"

#### Add your code to install here.
mkdir -p "${TARGET_DIR}"

# These are set to non-readable
# Causes installer to fail when removing temp dir
chmod -R a+rw ${EXTRACTED_DIR_CALLED}/syn/*

cp -rL ${EXTRACTED_DIR_CALLED}/* "${TARGET_DIR}"
cd "${TARGET_DIR}"

cp ../../signalp.patch ./
patch signalp signalp.patch
sed -i "s~/usr/opt/www/pub/CBS/services/SignalP-4.1/signalp-4.1~${WORK_DIR}/${TARGET_DIR}~" ./signalp

chmod -R a+r ./*
chmod a+rx signalp bin/* lib/*

#nb we delete WORKDIR using a trap command in register-base.sh
chmod -R u+w ${WORK_DIR}/${EXTRACTED_DIR_CALLED}
rm -r ${WORK_DIR}/${EXTRACTED_DIR_CALLED}

# Don't change the next two lines
echo "Finished registering ${EXTRACTED_DIR_CALLED}."
echo "Testing installation..."

# Add a command that uses the test dataset included with the package.
# both TEST_RESULT and TEST_RETCODE should be set.
# If you REALLY have to skip tests, set TEST_RETCODE=0.

export PATH=${WORK_DIR}/${TARGET_DIR}:$PATH
# If you get a non-zero exit code, the test will fail.
TEST_RESULT=$(signalp "${WORK_DIR}/${TARGET_DIR}/test/euk10.fsa")
TEST_RETCODE=$?

if [ "$TEST_RETCODE" -eq 0 ]; then
    echo "SignalP test successful"
else
    echo "SignalP test failed with exit code $TEST_RETCODE"
fi

echo "To include SignalP in your Nextflow run, remember to add:"
echo "--signalp_path \"${WORK_DIR}/${TARGET_DIR}/\""
