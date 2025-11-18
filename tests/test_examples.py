import sys
from pathlib import Path

# Add parent directory to path
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

from examples import boat_reactor_example, coated_wall_reactor_example


def test_boat_reactor_example():
    boat_reactor_example.main()


def test_coated_wall_reactor_example():
    coated_wall_reactor_example.main()
