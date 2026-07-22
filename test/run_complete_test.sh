#!/bin/bash
# 完整的一键测试脚本
# 用法: bash run_complete_test.sh

set -e  # 遇到错误立即退出

echo "=========================================="
echo " dockingML 完整流程测试"
echo "=========================================="

# 检查GROMACS
echo ""
echo "[1] 检查GROMACS..."
if ! command -v gmx &> /dev/null; then
    echo "错误: 未找到GROMACS (gmx命令)"
    echo "请安装GROMACS: sudo apt-get install gromacs"
    exit 1
fi
echo "✓ GROMACS已安装: $(which gmx)"

# 创建测试目录
TEST_DIR="test/test_output_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"
echo ""
echo "[2] 工作目录: $(pwd)"

# 复制测试PDB
echo ""
echo "[3] 准备测试PDB..."
cp ../test_data/test_protein.pdb input.pdb
echo "✓ 已复制测试PDB文件"

# 步骤1: pdb2gmx
echo ""
echo "[4] 运行 pdb2gmx..."
gmx pdb2gmx -f input.pdb -o processed.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh << EOF
1
EOF
if [ $? -eq 0 ]; then
    echo "✓ pdb2gmx 成功"
else
    echo "✗ pdb2gmx 失败"
    exit 1
fi

# 步骤2: editconf
echo ""
echo "[5] 运行 editconf (添加盒子)..."
gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic
if [ $? -eq 0 ]; then
    echo "✓ editconf 成功"
else
    echo "✗ editconf 失败"
    exit 1
fi

# 步骤3: solvate (加溶剂) - 关键步骤
echo ""
echo "[6] 运行 solvate (加溶剂)..."
echo "  这是之前出错的地方，仔细检查..."

# 查找spc903.gro
SPC_FILE=""
for path in "../automd/data/spc903.gro" "/usr/share/gromacs/top/spc903.gro" "/usr/local/share/gromacs/top/spc903.gro"; do
    if [ -f "$path" ]; then
        SPC_FILE="$path"
        echo "  找到spc文件: $SPC_FILE"
        break
    fi
done

if [ -z "$SPC_FILE" ]; then
    echo "  未找到spc903.gro文件，使用gmx默认溶剂盒子..."
    gmx solvate -cp boxed.gro -o solvated.gro -p topol.top
else
    gmx solvate -cp boxed.gro -cs "$SPC_FILE" -o solvated.gro -p topol.top
fi

if [ $? -eq 0 ]; then
    echo "✓ solvate 成功"
else
    echo "✗ solvate 失败"
    echo ""
    echo " [诊断] 请检查:"
    echo "   1. spc903.gro文件是否存在"
    echo "   2. boxed.gro文件是否有效"
    echo "   3. topol.top文件是否有效"
    echo "   4. 磁盘空间是否足够"
    exit 1
fi

# 步骤4: genion (加离子)
echo ""
echo "[7] 运行 genion (加离子)..."
gmx grompp -f ../data/em_sol.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 100
echo "SOL" | gmx genion -s ions.tpr -p topol.top -o ionized.gro -neutral -conc 0.15
if [ $? -eq 0 ]; then
    echo "✓ genion 成功"
else
    echo "✗ genion 失败"
    exit 1
fi

# 步骤5: 能量最小化
echo ""
echo "[8] 运行能量最小化..."
gmx grompp -f ../data/em_sol.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 100
gmx mdrun -deffnm em -nt 1 -v
if [ $? -eq 0 ]; then
    echo "✓ 能量最小化成功"
else
    echo "✗ 能量最小化失败"
    exit 1
fi

echo ""
echo "=========================================="
echo " 测试完成！所有步骤都成功。"
echo "=========================================="
echo ""
echo "输出文件在: $(pwd)"
ls -lh

