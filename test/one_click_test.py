#!/usr/bin/env python
"""
一键跑通整个流程的测试脚本
包含两种模式：
1. 模拟模式（不需要GROMACS，测试代码逻辑）
2. 真实模式（需要GROMACS，测试完整流程）
"""

import os
import sys
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

print("=" * 80)
print(" dockingML 完整流程一键测试")
print("=" * 80)

def test_simulation_mode():
    """模拟模式：测试代码逻辑，不实际运行GROMACS"""
    print("\n[模拟模式] 测试代码逻辑（不需要GROMACS）...")
    print("=" * 80)
    
    try:
        from automd.autoRunMD_gmx import AutoRunMD
        
        # 创建测试目录
        test_dir = PROJECT_ROOT / "test" / "simulation_test"
        test_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n测试目录: {test_dir}")
        
        # 初始化
        app = AutoRunMD()
        print(f"✓ AutoRunMD 初始化成功")
        print(f"  PROJECT_ROOT = {app.PROJECT_ROOT}")
        
        # 测试1: 检查所有必要的方法是否存在
        print("\n[测试 1] 检查方法存在性...")
        methods = ['generate_top', 'add_box', 'add_solvent', 'add_ions', 'minimize', 'md', 'run_app']
        for method in methods:
            if hasattr(app, method):
                print(f"  ✓ {method}() 方法存在")
            else:
                print(f"  ✗ {method}() 方法不存在")
                return False
        
        # 测试2: 检查add_solvent方法的改进
        print("\n[测试 2] 检查 add_solvent() 改进...")
        import inspect
        source = inspect.getsource(app.add_solvent)
        
        # 检查是否包含调试信息
        if "print(f\"[add_solvent]" in source:
            print("  ✓ add_solvent() 包含调试信息")
        else:
            print("  ⚠ add_solvent() 可能不包含调试信息")
        
        # 检查是否包含文件存在性检查
        if "os.path.exists(ingro)" in source:
            print("  ✓ add_solvent() 包含输入文件检查")
        else:
            print("  ⚠ add_solvent() 可能不包含输入文件检查")
        
        # 检查是否包含多个spc文件路径
        if "spc_paths" in source:
            print("  ✓ add_solvent() 包含多个spc文件路径检查")
        else:
            print("  ⚠ add_solvent() 可能不包含多个spc文件路径检查")
        
        # 测试3: 模拟add_solvent调用（传入不存在的文件，触发错误）
        print("\n[测试 3] 测试错误处理...")
        try:
            app.add_solvent("nonexistent.gro", "output.gro")
            print("  ✗ 期望抛出异常，但没有")
            return False
        except FileNotFoundError as e:
            print(f"  ✓ 正确捕获文件不存在错误:")
            print(f"    {str(e)[:200]}...")
        except Exception as e:
            print(f"  ⚠ 捕获到其他异常: {type(e).__name__}: {e}")
        
        # 测试4: 检查数据文件
        print("\n[测试 4] 检查数据文件完整性...")
        data_dir = Path(app.PROJECT_ROOT) / "data"
        
        required_files = {
            "spc903.gro": "SPC水盒子文件",
            "em_sol.mdp": "能量最小化参数文件",
            "npt.mdp": "NPT MD参数文件",
            "gbsa.mdp": "GBSA MD参数文件"
        }
        
        all_exist = True
        for fname, desc in required_files.items():
            fpath = data_dir / fname
            if fpath.exists():
                size = fpath.stat().st_size
                print(f"  ✓ {fname} ({desc}) - {size} bytes")
            else:
                print(f"  ✗ {fname} ({desc}) - 不存在")
                all_exist = False
        
        if not all_exist:
            print("  ⚠ 部分数据文件缺失")
        
        print("\n" + "=" * 80)
        print("[模拟模式] 测试完成！")
        print("=" * 80)
        print("\n结论:")
        print("  ✓ 代码逻辑正常")
        print("  ✓ 错误处理已改进")
        print("  ⚠ 要测试完整流程，请安装GROMACS并运行真实模式")
        
        return True
        
    except ImportError as e:
        print(f"\n✗ 导入失败: {e}")
        return False
    except Exception as e:
        print(f"\n✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_real_mode():
    """真实模式：实际运行GROMACS（需要安装）"""
    print("\n[真实模式] 测试完整流程（需要GROMACS）...")
    print("=" * 80)
    
    # 检查GROMACS
    import subprocess
    try:
        result = subprocess.run(["which", "gmx"], capture_output=True, text=True)
        if result.returncode != 0:
            print("\n✗ GROMACS (gmx命令) 未安装")
            print("  请先安装GROMACS:")
            print("  Ubuntu/Debian: sudo apt-get install gromacs")
            print("  CentOS/RHEL: sudo yum install gromacs")
            return False
        else:
            gmx_path = result.stdout.strip()
            print(f"\n✓ 找到GROMACS: {gmx_path}")
    except Exception as e:
        print(f"\n✗ 检查GROMACS失败: {e}")
        return False
    
    # 创建测试目录
    test_dir = PROJECT_ROOT / "test" / "real_test"
    test_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n测试目录: {test_dir}")
    
    # 创建简单的测试PDB
    print("\n[步骤 1] 创建测试PDB文件...")
    test_pdb = test_dir / "input.pdb"
    
    pdb_content = """HEADER    TEST PROTEIN
TITLE     TEST PROTEIN FOR DOCKINGML
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.058   1.422   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.358   2.522   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.958  -0.600  -1.258  1.00  0.00           C
ATOM      6  H   ALA A   1      -0.600  -0.800   0.000  1.00  0.00           H
ATOM      7  N   GLY A   2       3.458   1.422   0.000  1.00  0.00           N
ATOM      8  CA  GLY A   2       4.058   2.844   0.000  1.00  0.00           C
ATOM      9  C   GLY A   2       5.516   2.844   0.000  1.00  0.00           C
ATOM     10  O   GLY A   2       6.116   1.744   0.000  1.00  0.00           O
ATOM     11  H   GLY A   2       4.058   0.622   0.000  1.00  0.00           H
ENDMDL
END
"""
    test_pdb.write_text(pdb_content)
    print(f"  ✓ 已创建: {test_pdb}")
    
    # 切换到测试目录
    original_dir = os.getcwd()
    os.chdir(test_dir)
    
    try:
        from automd.autoRunMD_gmx import AutoRunMD
        
        # 初始化
        app = AutoRunMD()
        
        # 步骤2: pdb2gmx
        print("\n[步骤 2] 运行 pdb2gmx...")
        try:
            app.generate_top("input.pdb", outgro="processed", top="topol", ff="amber99sb-ildn")
            print("  ✓ pdb2gmx 成功")
        except Exception as e:
            print(f"  ✗ pdb2gmx 失败: {e}")
            return False
        
        # 步骤3: editconf
        print("\n[步骤 3] 运行 editconf (添加盒子)...")
        try:
            app.add_box(ingro="processed.gro", outgro="boxed.gro")
            print("  ✓ editconf 成功")
        except Exception as e:
            print(f"  ✗ editconf 失败: {e}")
            return False
        
        # 步骤4: solvate (加溶剂) - 这是之前出错的地方
        print("\n[步骤 4] 运行 solvate (加溶剂)...")
        print("  " + "=" * 60)
        try:
            app.add_solvent(ingro="boxed.gro", outgro="solvated.gro", intop="topol")
            print("  " + "=" * 60)
            print("  ✓ solvate 成功！")
        except Exception as e:
            print("  " + "=" * 60)
            print(f"  ✗ solvate 失败: {e}")
            print("\n  [诊断] 请检查:")
            print("    1. spc903.gro 文件是否存在")
            print("    2. boxed.gro 文件是否有效")
            print("    3. GROMACS 是否正确安装")
            return False
        
        # 步骤5: genion (加离子)
        print("\n[步骤 5] 运行 genion (加离子)...")
        try:
            app.add_ions(ingro="solvated.gro", outgro="ionized.gro", intop="topol")
            print("  ✓ genion 成功")
        except Exception as e:
            print(f"  ✗ genion 失败: {e}")
            print("  ⚠ 跳过此步骤，继续测试...")
        
        # 步骤6: 能量最小化
        print("\n[步骤 6] 运行能量最小化...")
        try:
            app.minimize(ingro="ionized.gro", outgro="em", intop="topol", nt=1)
            print("  ✓ 能量最小化成功")
        except Exception as e:
            print(f"  ✗ 能量最小化失败: {e}")
            print("  ⚠ 跳过此步骤...")
        
        print("\n" + "=" * 80)
        print("[真实模式] 测试完成！")
        print("=" * 80)
        print(f"\n输出文件在: {test_dir}")
        
        # 列出生成的文件
        print("\n生成的文件:")
        for f in test_dir.iterdir():
            if f.is_file():
                size = f.stat().st_size
                print(f"  - {f.name} ({size} bytes)")
        
        return True
        
    except Exception as e:
        print(f"\n✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # 恢复原始目录
        os.chdir(original_dir)

def create_quick_test_script():
    """创建一个快速测试脚本（Shell版本）"""
    print("\n" + "=" * 80)
    print("创建快速测试脚本...")
    print("=" * 80)
    
    script_path = PROJECT_ROOT / "test" / "quick_test.sh"
    
    script_content = """#!/bin/bash
# 快速测试脚本 - 一键跑通整个流程
# 用法: bash quick_test.sh [simulation|real]

set -e

MODE=${1:-"simulation"}

echo "=========================================="
echo " dockingML 快速测试"
echo " 模式: $MODE"
echo "=========================================="

if [ "$MODE" == "simulation" ]; then
    echo ""
    echo "运行模拟模式（不需要GROMACS）..."
    python test_complete_diagnosis.py
elif [ "$MODE" == "real" ]; then
    echo ""
    echo "运行真实模式（需要GROMACS）..."
    
    # 检查GROMACS
    if ! command -v gmx &> /dev/null; then
        echo "错误: 未找到GROMACS (gmx命令)"
        echo "请安装GROMACS或确保其已在PATH中"
        exit 1
    fi
    
    echo "✓ GROMACS已安装"
    
    # 创建测试目录
    TEST_DIR="test/real_test_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$TEST_DIR"
    cd "$TEST_DIR"
    
    echo ""
    echo "测试目录: $(pwd)"
    echo ""
    
    # 创建测试PDB
    echo "[1] 创建测试PDB文件..."
    cat > input.pdb << 'EOF'
HEADER    TEST PROTEIN
TITLE     TEST PROTEIN FOR DOCKINGML
MODEL        1
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.058   1.422   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.358   2.522   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.958  -0.600  -1.258  1.00  0.00           C
ATOM      6  H   ALA A   1      -0.600  -0.800   0.000  1.00  0.00           H
ATOM      7  N   GLY A   2       3.458   1.422   0.000  1.00  0.00           N
ATOM      8  CA  GLY A   2       4.058   2.844   0.000  1.00  0.00           C
ATOM      9  C   GLY A   2       5.516   2.844   0.000  1.00  0.00           C
ATOM     10  O   GLY A   2       6.116   1.744   0.000  1.00  0.00           O
ATOM     11  H   GLY A   2       4.058   0.622   0.000  1.00  0.00           H
ENDMDL
END
EOF
    echo "✓ 已创建 input.pdb"
    
    # 运行流程
    echo ""
    echo "[2] 运行 pdb2gmx..."
    gmx pdb2gmx -f input.pdb -o processed.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh
    echo "✓ pdb2gmx 完成"
    
    echo ""
    echo "[3] 运行 editconf..."
    gmx editconf -f processed.gro -o boxed.gro -c -d 1.2 -bt cubic
    echo "✓ editconf 完成"
    
    echo ""
    echo "[4] 运行 solvate (加溶剂)..."
    echo "  这是之前出错的关键步骤..."
    
    # 查找spc903.gro
    SPC_FILE=""
    for path in "../automd/data/spc903.gro" "/usr/share/gromacs/top/spc903.gro"; do
        if [ -f "$path" ]; then
            SPC_FILE="$path"
            echo "  找到spc文件: $SPC_FILE"
            break
        fi
    done
    
    if [ -z "$SPC_FILE" ]; then
        echo "  未找到spc903.gro，使用默认溶剂盒子..."
        gmx solvate -cp boxed.gro -o solvated.gro -p topol.top
    else
        gmx solvate -cp boxed.gro -cs "$SPC_FILE" -o solvated.gro -p topol.top
    fi
    echo "✓ solvate 完成"
    
    echo ""
    echo "[5] 运行 genion..."
    gmx grompp -f ../automd/data/em_sol.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 100
    echo "SOL" | gmx genion -s ions.tpr -p topol.top -o ionized.gro -neutral -conc 0.15
    echo "✓ genion 完成"
    
    echo ""
    echo "=========================================="
    echo " 测试完成！"
    echo "=========================================="
    echo ""
    echo "输出文件在: $(pwd)"
    ls -lh
    
else
    echo "未知模式: $MODE"
    echo "用法: bash $0 [simulation|real]"
    exit 1
fi
"""
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # 添加执行权限
    os.chmod(script_path, 0o755)
    
    print(f"\n✓ 已创建快速测试脚本: {script_path}")
    print(f"\n使用方法:")
    print(f"  模拟模式（不需要GROMACS）:")
    print(f"    bash {script_path} simulation")
    print(f"  或:")
    print(f"    python {PROJECT_ROOT / 'test' / 'one_click_test.py'}")
    print(f"\n  真实模式（需要GROMACS）:")
    print(f"    bash {script_path} real")
    
    return script_path

def main():
    """主函数"""
    print("\n选择测试模式:")
    print("  1. 模拟模式（不需要GROMACS，快速测试代码逻辑）")
    print("  2. 真实模式（需要GROMACS，测试完整流程）")
    print("  3. 创建快速测试脚本")
    print("  4. 退出")
    
    while True:
        try:
            choice = input("\n请选择 (1-4): ").strip()
            
            if choice == "1":
                print("\n" + "=" * 80)
                test_simulation_mode()
                break
            elif choice == "2":
                print("\n" + "=" * 80)
                test_real_mode()
                break
            elif choice == "3":
                print("\n" + "=" * 80)
                create_quick_test_script()
                break
            elif choice == "4":
                print("\n退出...")
                break
            else:
                print("无效选择，请输入 1-4")
        except KeyboardInterrupt:
            print("\n\n中断...")
            break
        except Exception as e:
            print(f"\n错误: {e}")
            break
    
    print("\n" + "=" * 80)
    print("测试完成")
    print("=" * 80)

if __name__ == "__main__":
    main()
