use std::ops::MulAssign;

use num_traits::One;

use super::field_element::FieldElement;
use super::poly::Polynomial;
use super::traits::FiniteField;

#[derive(Debug, Clone, PartialEq)]
pub struct Leaf<FF: FiniteField + MulAssign<FieldElement> + 'static> {
    pub(crate) points: Vec<FF>,
    zerofier: Polynomial<'static, FF>,
}

impl<FF> Leaf<FF>
where
    FF: FiniteField + MulAssign<FieldElement> + 'static,
{
    pub fn new(points: Vec<FF>) -> Self {
        let zerofier = Polynomial::zerofier(&points);
        Self { points, zerofier }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Branch<FF: FiniteField + MulAssign<FieldElement> + 'static> {
    zerofier: Polynomial<'static, FF>,
    pub(crate) left: ZerofierTree<FF>,
    pub(crate) right: ZerofierTree<FF>,
}

impl<FF> Branch<FF>
where
    FF: FiniteField + MulAssign<FieldElement> + 'static,
{
    pub fn new(left: ZerofierTree<FF>, right: ZerofierTree<FF>) -> Self {
        let zerofier = left.zerofier().multiply(&right.zerofier());
        Self {
            zerofier,
            left,
            right,
        }
    }
}

/// A zerofier tree is a balanced binary tree of vanishing polynomials.
/// Conceptually, every leaf corresponds to a single point, and the value of
/// that leaf is the monic linear polynomial that evaluates to zero there and
/// no-where else. Every non-leaf node is the product of its two children.
/// In practice, it makes sense to truncate the tree depth, in which case every
/// leaf contains a chunk of points whose size is upper-bounded and more or less
/// equal to some constant threshold.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum ZerofierTree<FF: FiniteField + MulAssign<FieldElement> + 'static> {
    Leaf(Leaf<FF>),
    Branch(Box<Branch<FF>>),
    #[default]
    Padding,
}

impl<FF: FiniteField + MulAssign<FieldElement> + 'static> ZerofierTree<FF> {
    /// Regulates the depth at which the tree is truncated. Phrased differently,
    /// regulates the number of points contained by each leaf.
    const RECURSION_CUTOFF_THRESHOLD: usize = 16;

    pub fn new_from_domain(domain: &[FF]) -> Self {
        Self::build(domain)
    }

    fn build(domain: &[FF]) -> Self {
        if domain.is_empty() {
            return ZerofierTree::Padding;
        }
        if domain.len() <= Self::RECURSION_CUTOFF_THRESHOLD {
            return ZerofierTree::Leaf(Leaf::new(domain.to_vec()));
        }

        let mid = domain.len() / 2;
        let left = Self::build(&domain[..mid]);
        let right = Self::build(&domain[mid..]);
        Self::combine(left, right)
    }

    fn combine(left: Self, right: Self) -> Self {
        match (left, right) {
            (ZerofierTree::Padding, r) => r,
            (l, ZerofierTree::Padding) => l,
            (l, r) => ZerofierTree::Branch(Box::new(Branch::new(l, r))),
        }
    }

    pub fn zerofier(&self) -> Polynomial<'static, FF> {
        match self {
            ZerofierTree::Leaf(leaf) => leaf.zerofier.clone(),
            ZerofierTree::Branch(branch) => branch.zerofier.clone(),
            ZerofierTree::Padding => Polynomial::one(),
        }
    }
}

#[cfg(test)]
mod test {
    use num_traits::ConstZero;
    use num_traits::Zero;
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest_arbitrary_interop::arb;
    use test_strategy::proptest;

    use crate::prelude::FieldElement;
    use crate::prelude::Polynomial;
    use crate::zerofier_tree::ZerofierTree;

    fn domain_points(range: std::ops::Range<usize>) -> Vec<FieldElement> {
        range
            .map(|i| FieldElement::from(i as u32))
            .collect::<Vec<_>>()
    }

    #[test]
    fn zerofier_tree_can_be_empty() {
        ZerofierTree::<FieldElement>::new_from_domain(&[]);
    }

    #[test]
    fn new_from_domain_returns_leaf() {
        let count = ZerofierTree::<FieldElement>::RECURSION_CUTOFF_THRESHOLD;
        let points = domain_points(0..count);
        match ZerofierTree::new_from_domain(&points) {
            ZerofierTree::Leaf(leaf) => {
                assert_eq!(leaf.points, points);
                assert_eq!(Polynomial::zerofier(&leaf.points), leaf.zerofier);
            }
            _ => panic!("expected leaf"),
        }
    }

    #[test]
    fn new_from_domain_returns_branch() {
        let count =
            ZerofierTree::<FieldElement>::RECURSION_CUTOFF_THRESHOLD + 1;
        let points = domain_points(0..count);
        match ZerofierTree::new_from_domain(&points) {
            ZerofierTree::Branch(branch) => {
                assert_ne!(branch.left, ZerofierTree::Padding);
                assert_ne!(branch.right, ZerofierTree::Padding);
                assert_eq!(
                    branch.zerofier,
                    branch.left.zerofier().multiply(&branch.right.zerofier())
                );
            }
            _ => panic!("expected branch"),
        }
    }
    #[proptest]
    fn zerofier_tree_root_is_multiple_of_children(
        #[strategy(vec(arb(), 2*ZerofierTree::<FieldElement>::RECURSION_CUTOFF_THRESHOLD))]
        points: Vec<FieldElement>,
    ) {
        let zerofier_tree = ZerofierTree::new_from_domain(&points);
        let ZerofierTree::Branch(ref branch) = &zerofier_tree else {
            panic!("not enough leafs");
        };
        prop_assert_eq!(
            Polynomial::zero(),
            zerofier_tree.zerofier().reduce(&branch.left.zerofier())
        );
        prop_assert_eq!(
            Polynomial::zero(),
            zerofier_tree.zerofier().reduce(&branch.right.zerofier())
        );
    }

    #[proptest]
    fn zerofier_tree_root_has_right_degree(
        #[strategy(vec(arb(), 1..(1<<10)))] points: Vec<FieldElement>,
    ) {
        let zerofier_tree = ZerofierTree::new_from_domain(&points);
        prop_assert_eq!(
            points.len(),
            zerofier_tree.zerofier().degree() as usize
        );
    }

    #[proptest]
    fn zerofier_tree_root_zerofies(
        #[strategy(vec(arb(), 1..(1<<10)))] points: Vec<FieldElement>,
        #[strategy(0usize..#points.len())] index: usize,
    ) {
        let zerofier_tree = ZerofierTree::new_from_domain(&points);
        prop_assert_eq!(
            FieldElement::ZERO,
            zerofier_tree.zerofier().evaluate(points[index])
        );
    }

    #[proptest]
    fn zerofier_tree_and_polynomial_agree_on_zerofiers(
        #[strategy(vec(arb(), 1..(1<<10)))] points: Vec<FieldElement>,
    ) {
        let zerofier_tree = ZerofierTree::new_from_domain(&points);
        let polynomial_zerofier = Polynomial::zerofier(&points);
        prop_assert_eq!(polynomial_zerofier, zerofier_tree.zerofier());
    }
}
