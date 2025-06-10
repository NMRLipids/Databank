import os
import sys
import subprocess



def main():
    source_dir = os.path.abspath(os.path.dirname(__file__))
    # Project root is three levels up
    repo_root = os.path.abspath(os.path.join(source_dir, '..', '..', '..'))

    scripts_dir = os.path.relpath(os.path.join(repo_root, 'Scripts'), source_dir)
    auto_dir = 'auto_gen'
    template_dir = '_templates/apidoc'
    excluded_patterns = ['*tests*','*init*']


    cmd = [
        'sphinx-apidoc',     
        '-o', auto_dir,       # relative output directory
        '--separate',         # separate .rst per module
        '--implicit-namespaces',  # support namespace packages
        '-d', '1',            # max depth = 1
        '--templatedir', template_dir,
        '--tocfile','Scripts',   # relative template directory
        scripts_dir, *excluded_patterns           
    ]

    # Display info and execute within source_dir
    print('Running command:')
    print(' '.join(cmd))

    try:
        subprocess.run(cmd, cwd=source_dir)
        print('sphinx-apidoc completed successfully')
    except subprocess.CalledProcessError as e:
        print(f'sphinx-apidoc failed with exit code {e.returncode}', file=sys.stderr)
        sys.exit(e.returncode)




if __name__ == '__main__':
    main()