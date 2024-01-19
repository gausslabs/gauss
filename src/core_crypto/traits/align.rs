//! We cannot impl the get_col_iter for Matrix at the trait level since it is not object-safe. 
//! However, we can impl this method for common 2d-array like data strcture without much effort.
//! 


use aligned_vec::AVec;

pub trait Matrix: AsRef<[Self::R]>
{
    type  MatElement;
    type R:Row<Element = Self::MatElement>;

    fn get_row(&self,row_idx:usize)->&Self::R{
        &self.as_ref()[row_idx]
    }
    fn get_col_iter(&self,column_idx:usize)->impl Iterator<Item = &Self::MatElement>{
        self.as_ref().iter().map( move |r|{
            &r.as_ref()[column_idx]
        })
    }
    fn get_element(&self,row_idx:usize,column_idx:usize)->&Self::MatElement{
        &self.as_ref()[row_idx].as_ref()[column_idx]
    }
    
}



pub trait MatrixMut:Matrix+AsMut<[<Self as  Matrix>::R]>
where <Self as Matrix>::R:RowMut {
    fn get_row_mut(&mut self,row_idx:usize)->&mut [<Self as Matrix>::MatElement]{
        self.as_mut()[row_idx].as_mut()
    }

    #[inline]
    fn set(&mut self, row_idx:usize,column_idx:usize,val:<Self as Matrix>::MatElement){
       
       self.as_mut()[row_idx].as_mut()[column_idx]=val;

    }
    fn get_col_iter_mut(&mut self,column_idx:usize)->impl Iterator<Item = &mut Self::MatElement> {
        let res = self.as_mut().iter_mut().map(move |r| {
           r.as_mut().get_mut(column_idx).unwrap()
        });
           res
       }
       fn get_element_mut(&mut self,row_idx:usize,column_idx:usize)->&mut <Self as Matrix>::MatElement{
        &mut self.as_mut()[row_idx].as_mut()[column_idx]
    }
}

  

 ///This part is borrowed from zama fhe's containter type without modification
 ///https://github.com/zama-ai/tfhe-rs/blob/main/tfhe/src/core_crypto/commons/traits/container.rs
 /// 
pub trait Row :AsRef<[Self::Element]>{
    type Element;
}

pub trait RowMut: Row + AsMut<[<Self as Row> ::Element]> {}
impl<T>  Row for AVec<T> {
    type Element = T;
}
impl<T> Row for Vec<T> {
    type Element = T;
}

impl<T> RowMut for Vec<T> {}

impl<T> Row for &[T] {
    type Element = T;
}

impl<T> Row for &mut [T] {
    type Element = T;
}

impl<T> RowMut for &mut [T] {}

impl<T>  Row for aligned_vec::ABox<[T]> {
    type Element = T;
}

impl<T> RowMut for aligned_vec::ABox<[T]> {}

impl<T>  Row for Box<[T]> {
    type Element = T;
}

impl<T> RowMut for Box<[T]> {}
impl<T> RowMut for aligned_vec::AVec<T> {}


pub trait IntoRowOwned: Row + AsMut<[Self::Element]> {
    fn collect<I: Iterator<Item = Self::Element>>(iter: I) -> Self;
}

impl<T> IntoRowOwned for aligned_vec::ABox<[T]> {
    fn collect<I: Iterator<Item = Self::Element>>(iter: I) -> Self {
        aligned_vec::AVec::<T, _>::from_iter(0, iter).into_boxed_slice()
    }
}

///
/// 
/// 

type InnerMat<T> = AVec<AVec<T>>;
impl<T> Matrix for InnerMat<T>{
    type  MatElement=T;
    type R=AVec<T>;
}



impl<T> MatrixMut for InnerMat<T> {}

impl<T> Matrix for Vec<Vec<T>>{
    type  MatElement=T;
    type R=Vec<T>;
}
impl<T> MatrixMut for Vec<Vec<T>> {}

#[cfg(test)]
mod test{
   
    use super::*;
    use aligned_vec::avec;
    use crate::core_crypto::traits::align::MatrixMut;
    
    #[test]
    fn test(){
        let  v1 = avec![1_u64,2_u64];
        let  v2 = avec![3_u64,4_u64];
        let  mut v3 = avec![v1,v2];
        let  mut r = v3.get_row_mut(0);
        r[0]=2;
        println!("{:?}",v3);
        let v5 :Vec<_>=v3.get_col_iter_mut(0usize).collect();
        for i in v5 {
            *i=0;
        }
        println!("{:?}",v3);
    }
}
