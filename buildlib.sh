#!/bin/bash

ROOT_FOLDER=$(pwd)
DFT_LIB="dft_lib"
BUILD_DIR="build/release"
OPTION_FLAG=""

go_to_root_dir(){
  local root_dir=$ROOT_FOLDER
  cd $root_dir
}

prepare_subfolder(){
  local folder=$1
  local build_type=$2

  mkdir -p $folder
  cd $folder

  echo $folder $build_typep
  cmake -DCMAKE_BUILD_TYPE=$build_type ../../
}

prepare_folders(){
  echo "======================================"
  echo "  Starting folder configuration"
  echo "======================================"
  local root_dir=$ROOT_FOLDER
  parent_dir="$(dirname "$root_dir")"

  mkdir -p build

  cd build
  prepare_subfolder "release" "Release"
  go_to_root_dir

  cd build
  prepare_subfolder "debug" "Debug"
  go_to_root_dir
}

default_build_lib(){
  local build_dir=$1
  local option_flag=$2

  echo cmake --build $build_dir $option_flag
  cmake --build $build_dir $option_flag
}

build_lib(){
  echo "Runninf default build..."
  echo $ROOT_FOLDER
  cd ../$DFT_LIB

  echo "Making Library"
  echo cmake --build $BUILD_DIR $OPTION_FLAG
  cmake --build $BUILD_DIR $OPTION_FLAG

  cp core/
  cd $ROOT_FOLDER
}

run_tests(){
  local build_dir=$1
  cd $build_dir
  make gtests
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
          local build_dir=$BUILD_DIR
          local option_flag="--clean-first"
          default_build_lib $build_dir $option_flag
          ;;
        debug)
          local build_dir="build/debug"
          local option_flag=$OPTION_FLAG
          default_build_lib $build_dir $option_flag
          ;;
	      lib)
          local build_dir="build/release"
          local option_flag=$OPTION_FLAG

          echo "[Step 1/1] Building Library..."
          default_build_lib $build_dir $option_flag
          ;;
        test)
          local build_dir="build/release"
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
main $exec_mode