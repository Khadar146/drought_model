#!/bin/bash
# Build script for Somaliland Drought Analysis System
# Handles NetCDF linking automatically

echo "🏗️  Building Somaliland Drought Analysis System..."

# Set NetCDF paths (adjust these paths if needed)
NETCDF_INCLUDE="/Users/kdc/miniconda3/envs/drought_env/include"
NETCDF_LIB="/Users/kdc/miniconda3/envs/drought_env/lib"

# Build flags
BUILD_FLAGS="-I${NETCDF_INCLUDE} -L${NETCDF_LIB} -lnetcdff -lnetcdf"

echo "📦 NetCDF Include: $NETCDF_INCLUDE"
echo "📦 NetCDF Library: $NETCDF_LIB"
echo "🔧 Build flags: $BUILD_FLAGS"
echo ""

# Build the project
echo "⚙️  Compiling..."
fpm build --flag "$BUILD_FLAGS"

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Build successful!"
    echo ""
    echo "🚀 Usage examples:"
    echo "   fpm run DroughtAnalysis --flag '$BUILD_FLAGS' -- historical"
    echo "   fpm run DroughtAnalysis --flag '$BUILD_FLAGS' -- future"
    echo "   fpm run DroughtAnalysis --flag '$BUILD_FLAGS' -- comparison"
    echo ""
    echo "💡 Tip: Create an alias for easier running:"
    echo "   alias drought='fpm run DroughtAnalysis --flag \"$BUILD_FLAGS\" --'"
    echo ""
else
    echo "❌ Build failed!"
    exit 1
fi
