#!/bin/bash

DFT_LIB="dft_lib"
ROOT_FOLDER="$(pwd)/$DFT_LIB"
DFT_TEST="gtests"
BUILD_DIR="cmake-build"
DEBUG_BUILD_DIR="cmake-build/debug"
RELEASE_BUILD_DIR="cmake-build/release"
OPTION_FLAG=""

go_to_root_dir(){
  local root_dir=$ROOT_FOLDER
  cd $root_dir
}

prepare_subfolder(){
  local folder=$1
  local build_type=$2

  mkdir -p "$folder"
  cd "$folder"

  echo "$folder $build_typep"
  cmake -DCMAKE_BUILD_TYPE=$build_type ../../
}

prepare_folders(){
  echo "======================================"
  echo "  Starting folder configuration"
  echo "======================================"
  cd "$ROOT_FOLDER" || exit

  mkdir -p $BUILD_DIR

  cd $BUILD_DIR
  prepare_subfolder "release" "Release"
  go_to_root_dir

  cd $BUILD_DIR
  prepare_subfolder "debug" "Debug"
  go_to_root_dir
}

default_build_lib(){
  local build_dir=$1
  local option_flag=$2

  echo cmake --build $build_dir $option_flag
  cmake --build $build_dir $option_flag
}

run_tests(){
  local build_dir=$1
  cd $build_dir
  make $DFT_TEST
  cd $ROOT_FOLDER
}

switch_case_menu(){
  echo ""
  echo "======================================"
  echo "  classicalDFT library                "
  echo "======================================"

  local mode=$1

  case "$mode" in
        clean)
          local build_dir=$DEBUG_BUILD_DIR
          local option_flag="--clean-first"
          default_build_lib $build_dir $option_flag
          ;;
        debug)
          local build_dir=$DEBUG_BUILD_DIR
          local option_flag=$OPTION_FLAG
          default_build_lib $build_dir $option_flag
          ;;
	      lib)
          local build_dir=$RELEASE_BUILD_DIR
          local option_flag=$OPTION_FLAG

          echo "[Step 1/1] Building Library..."
          default_build_lib $build_dir $option_flag
          ;;
        test)
          local build_dir=$DEBUG_BUILD_DIR
          local option_flag=$OPTION_FLAG

          echo "[Step 1/2] Making Library..."
          default_build_lib $build_dir $option_flag

          echo "[Step 2/2] Running tests..."
          run_tests $build_dir
          ;;
        *)
          echo "$mode: Not implemented yet"
          ;;
  esac
}

main(){
  local input_param=$1

  prepare_folders

  switch_case_menu $input_param
}

# TODO: Improve this documentation and script
# @exec_mode: The type of execution we want to run:
#     - clean
#     - lib
#     - debug
#     - test
exec_mode=$1
main "$exec_mode"