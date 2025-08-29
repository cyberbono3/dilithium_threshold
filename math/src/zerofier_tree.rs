use std::collections::VecDeque;
use std::ops::MulAssign;

use num_traits::One;

use super::field_element::FieldElement;
use super::polynomial::Polynomial;
use super::traits::FiniteField;

#[derive(Debug, Clone, PartialEq)]
pub struct Leaf<'c, FF: FiniteField + MulAssign<FieldElement>> {
    pub(crate) points: Vec<FF>,
    zerofier: Polynomial<'c, FF>,
}

impl<FF> Leaf<'static, FF>
where
    FF: FiniteField + MulAssign<FieldElement>,
{
    pub fn new(points: Vec<FF>) -> Leaf<'static, FF> {
        let zerofier = Polynomial::zerofier(&points);
        Self { points, zerofier }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Branch<'c, FF: FiniteField + MulAssign<FieldElement>> {
    zerofier: Polynomial<'c, FF>,
    pub(crate) left: ZerofierTree<'c, FF>,
    pub(crate) right: ZerofierTree<'c, FF>,
}

impl<'c, FF> Branch<'c, FF>
where
    FF: FiniteField + MulAssign<FieldElement> + 'static,
{
    pub fn new(
        left: ZerofierTree<'c, FF>,
        right: ZerofierTree<'c, FF>,
    ) -> Self {
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
#[derive(Debug, Clone, PartialEq)]
pub enum ZerofierTree<'c, FF: FiniteField + MulAssign<FieldElement>> {
    Leaf(Leaf<'c, FF>),
    Branch(Box<Branch<'c, FF>>),
    Padding,
}

impl<FF: FiniteField + MulAssign<FieldElement>> ZerofierTree<'static, FF> {
    /// Regulates the depth at which the tree is truncated. Phrased differently,
    /// regulates the number of points contained by each leaf.
    const RECURSION_CUTOFF_THRESHOLD: usize = 16;


    pub fn new_from_domain(domain: &[FF]) -> Self {
        // Build initial leaves (no padding leaves).
        let mut nodes: VecDeque<_> = domain
            .chunks(Self::RECURSION_CUTOFF_THRESHOLD)
            .map(|chunk| ZerofierTree::Leaf(Leaf::new(chunk.to_vec())))
            .collect();

        if nodes.is_empty() {
            return ZerofierTree::Padding;
        }

        // Pair level-by-level; carry lone nodes forward.
        while nodes.len() > 1 {
            let mut next = VecDeque::with_capacity((nodes.len() + 1) / 2);
            while let Some(left) = nodes.pop_front() {
                if let Some(right) = nodes.pop_front() {
                    match (left, right) {
                        (ZerofierTree::Padding, r) => next.push_back(r),
                        (l, ZerofierTree::Padding) => next.push_back(l),
                        (l, r) => {
                            let node = Branch::new(l, r);
                            next.push_back(ZerofierTree::Branch(Box::new(
                                node,
                            )));
                        }
                    }
                } else {
                    // Odd count: carry the last node forward unchanged.
                    next.push_back(left);
                }
            }
            nodes = next;
        }

        nodes.pop_front().unwrap()
    }
}

impl<'c, FF> ZerofierTree<'c, FF>
where
    FF: FiniteField + MulAssign<FieldElement> + 'static,
{
    pub fn zerofier(&self) -> Polynomial<'c, FF> {
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

    #[test]
    fn zerofier_tree_can_be_empty() {
        ZerofierTree::<FieldElement>::new_from_domain(&[]);
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
