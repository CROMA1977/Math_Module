#pragma once

using namespace std;

class CVectorReference2FLOAT;
class CVectorReference3FLOAT;
class CVectorReference4FLOAT;

class CVector2FLOAT;
class CVector3FLOAT;
class CVector4FLOAT;

class CMatrix2FLOAT;
class CMatrix3FLOAT;
class CMatrix4FLOAT;

/// 二維的向量型態 2x1
class  CVectorReference2FLOAT {
public:
			float						&m_x, &m_y;											///< 空間座標
										/// 建構式
										CVectorReference2FLOAT( float & Value0, float & Value1 );
										/// 建構式 Array[2] = {m_x, m_y}
										CVectorReference2FLOAT( float * Array );
										/// 建構式 Array[2] = {m_x, m_y}
										CVectorReference2FLOAT( CVector2FLOAT & Vector);
										/// 拷貝建構式
										CVectorReference2FLOAT( const CVectorReference2FLOAT & Vector );

										/// 解構式
	virtual							   ~CVectorReference2FLOAT();
			//////////////////////////////////////////////////////////////////////////
										/// 設定向量 
			CVectorReference2FLOAT	&	Set(const float x, const float y );
										/// 設定向量 Array[2] = {m_x, m_y}
			CVectorReference2FLOAT	&	Set(const float * Array );
										/// 設定向量 V = V0 + (V1-V0) * t 
			CVectorReference2FLOAT	&	Set(const CVectorReference2FLOAT & V0, const CVectorReference2FLOAT & V1, const float t );
										/// 計算兩向量的 Dot (Ax * Bx + Ay * By)
			float						Dot(const CVectorReference2FLOAT & Vector ) const;
										/// 計算兩向量的方向 逆時針 > 0 順時針 < 0
			float						Cross(const CVectorReference2FLOAT & Vector ) const;
										/// 使兩向量交換 
			void						Swap(CVectorReference2FLOAT & Vector );
										/// 取得兩向量的弳度
			float						GetRadian(const CVectorReference2FLOAT & Vector) const;
										/// 取得兩向量的角度
			float						GetAngle( const CVectorReference2FLOAT & Vector ) const;
										/// 取得向量的長度
			float						GetLength( void ) const;
										/// 取得向量長度的平方
			float						GetLengthSquare( void ) const;

			//////////////////////////////////////////////////////////////////////////
										/// 取得單位化的向量
			CVector2FLOAT 				GetNormalize( void ) const;
										/// 取得反轉的向量 V[2] = {-m_x, -m_y}
			CVector2FLOAT 				GetReverse( void ) const;
			//////////////////////////////////////////////////////////////////////////
										/// 使向量單位化
			CVectorReference2FLOAT	&	Normalize( void );
										/// 使向量反轉 
			CVectorReference2FLOAT	&	Reverse( void );
										/// 檢查向量是否為零
			bool						IsZero( void ) const;
			////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 + 運算子 V[2] = [m_x + Vector.m_x, m_y + Vector.m_y]
			CVector2FLOAT				operator+( const CVectorReference2FLOAT & Vector ) const;
										/// 覆載對向量的 - 運算子 V[2] = [m_x - Vector.m_x, m_y - Vector.m_y]
			CVector2FLOAT				operator-( const CVectorReference2FLOAT & Vector ) const;
										/// 覆載對實數的 * 運算子 V[2] = [m_x * Value, m_y * Value]
			CVector2FLOAT				operator*( const float Value ) const;
										/// 覆載對實數的 / 運算子 V[2] = [m_x / Value, m_y / Value] 
			CVector2FLOAT				operator/( const float Value ) const;

										/// 覆載 = 運算子 V = Vector
			CVectorReference2FLOAT	&	operator=( const CVectorReference2FLOAT & Vector );
										/// 覆載 = 運算子 V = Vector
			CVectorReference2FLOAT	&	operator=(const CVector2FLOAT& Vector);

										/// 覆載對向量的 += 運算子 [m_x, m_y] = [m_x + Vector.m_x, m_y + Vector.m_y]
			CVectorReference2FLOAT	&	operator+=( const CVectorReference2FLOAT & Vector );
										/// 覆載對向量的 -= 運算子 [m_x, m_y] = [m_x - Vector.m_x, m_y - Vector.m_y] 
			CVectorReference2FLOAT	&	operator-=( const CVectorReference2FLOAT & Vector );
										/// 覆載對實數的 *= 運算子 [m_x, m_y] = [m_x * Value, m_y * Value]
			CVectorReference2FLOAT	&	operator*=( const float Value );
										/// 覆載對實數的 /= 運算子 [m_x, m_y] = [m_x / Value, m_y / Value]
			CVectorReference2FLOAT	&	operator/=( const float Value );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 == 運算子 dx * dx + dy * dy <  0.000001
			bool						operator==( const CVectorReference2FLOAT & Vector ) const ;
										/// 覆載對向量的 != 運算子 dx * dx + dy * dy >= 0.000001
			bool						operator!=( const CVectorReference2FLOAT & Vector ) const ;
};

/// 三維的向量 3x1
class CVectorReference3FLOAT {
public:
			float		 				&m_x, &m_y, &m_z;					///< 在空間座標時使用的名稱

										/// 建構式
										CVectorReference3FLOAT( float & Value0, float & Value1, float & Value2 );
										/// 建構式 Array[3] = {m_x, m_y, m_z};
										CVectorReference3FLOAT( float * Array );
										/// 建構式 Array[3] = {m_x, m_y, m_z};
										CVectorReference3FLOAT( CVector3FLOAT & Vector);
										/// 拷貝建構式 										 
										CVectorReference3FLOAT( const CVectorReference3FLOAT & Vector );

										/// 解構式
	virtual							   ~CVectorReference3FLOAT();

			//////////////////////////////////////////////////////////////////////////
										/// 設定向量 
			CVectorReference3FLOAT	&	Set( const float x, const float y, const float z );
										/// 設定向量 Array[3] = {m_x, m_y, m_z};
			CVectorReference3FLOAT	&	Set( const float * Array );
										/// 設定向量 V = V0 + (V1-V0) * t 
			CVectorReference3FLOAT	&	Set( const CVectorReference3FLOAT & V0, const CVectorReference3FLOAT & V1, const float t );
										/// 計算兩向量的 Dot (Ax * Bx + Ay * By + Az * Bz)
			float						Dot( const CVectorReference3FLOAT & Vector ) const;
										/// 計算兩向量的 Cross 
			CVector3FLOAT				Cross( const CVectorReference3FLOAT & Vector ) const;
										/// 使兩向量交換 
			void						Swap( CVectorReference3FLOAT & Vector );
										/// 取得兩向量的弳度
			float						GetRadian(const CVectorReference3FLOAT & Vector) const;
										/// 取得兩向量的角度
			float						GetAngle( const CVectorReference3FLOAT & Vector ) const;
										/// 取得向量的長度
			float						GetLength( void ) const;
										/// 取得向量長度的平方
			float						GetLengthSquare( void ) const;

			//////////////////////////////////////////////////////////////////////////
										/// 取得單位化的向量
			CVector3FLOAT 				GetNormalize( void ) const;
										/// 取得反轉的向量 V[3] = {-m_x, -m_y, -m_z}
			CVector3FLOAT 				GetReverse( void ) const;
			//////////////////////////////////////////////////////////////////////////
										/// 使向量單位化
			CVectorReference3FLOAT	&	Normalize( void );
										/// 使向量反轉 
			CVectorReference3FLOAT	&	Reverse( void );
										/// 檢查向量是否為零
			bool						IsZero( void ) const;
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 + 運算子 V[3] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z]
			CVector3FLOAT				operator+( const CVectorReference3FLOAT & Vector ) const;
										/// 覆載對向量的 - 運算子 V[3] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z]
			CVector3FLOAT				operator-( const CVectorReference3FLOAT & Vector ) const;
										/// 覆載對實數的 * 運算子 V[3] = [m_x * Value, m_y * Value, m_z * Value]
			CVector3FLOAT				operator*( const float Value ) const;
										/// 覆載對實數的 / 運算子 V[3] = [m_x / Value, m_y / Value, m_z / Value] 
			CVector3FLOAT				operator/( const float Value ) const;

										/// 覆載 = 運算子 V = Vector
			CVectorReference3FLOAT	&	operator=( const CVectorReference3FLOAT & Vector );
										/// 覆載 = 運算子 V = Vector
			CVectorReference3FLOAT	&	operator=( const CVector3FLOAT & Vector);

										/// 覆載對向量的 += 運算子 [m_x, m_y, m_z] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z]
			CVectorReference3FLOAT	&	operator+=( const CVectorReference3FLOAT & Vector );
										/// 覆載對向量的 -= 運算子 [m_x, m_y, m_z] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z] 
			CVectorReference3FLOAT	&	operator-=( const CVectorReference3FLOAT & Vector );

										/// 覆載對實數的 *= 運算子 [m_x, m_y, m_z] = [m_x * Value, m_y * Value, m_z * Value]
			CVectorReference3FLOAT	&	operator*=( const float Value );
										/// 覆載對實數的 /= 運算子 [m_x, m_y, m_z] = [m_x / Value, m_y / Value, m_z / Value] 
			CVectorReference3FLOAT	&	operator/=( const float Value );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 == 運算子 dx * dx + dy * dy + dz * dz <  0.000001
			bool						operator==( const CVectorReference3FLOAT  & Vector ) const;
										/// 覆載對向量的 != 運算子 dx * dx + dy * dy + dz * dz >= 0.000001
			bool						operator!=( const CVectorReference3FLOAT  & Vector ) const;

};

/// 四維的向量參考 4x1
class  CVectorReference4FLOAT {
public:
			float		 				&m_x, &m_y, &m_z, &m_w;		///< 在空間座標時使用的名稱

										/// 建構式
										CVectorReference4FLOAT( float & Value0, float & Value1, float & Value2, float & Value3 );
										/// 建構式 Array[4] = {m_x, m_y, m_z, m_w};
										CVectorReference4FLOAT( float * Array );
										/// 建構式 Array[4] = {m_x, m_y, m_z, m_w};
										CVectorReference4FLOAT( CVector4FLOAT & Vector);
										/// 拷貝建構式 								 
										CVectorReference4FLOAT( const CVectorReference4FLOAT & Vector );
										/// 解構式
	virtual							   ~CVectorReference4FLOAT();

			//////////////////////////////////////////////////////////////////////////
										/// 設定向量
			CVectorReference4FLOAT	&	Set( const float x, const float y, const float z, const float w );
										/// 設定向量 Array[4] = {m_x, m_y, m_z, m_w}
			CVectorReference4FLOAT	&	Set( const float * Array );									 
										/// 使兩向量交換 
			void						Swap( CVectorReference4FLOAT & Vector );		
			//////////////////////////////////////////////////////////////////////////
										/// 取得反轉的向量 V[4] = {-m_x, -m_y, -m_z, -m_w}
			CVector4FLOAT				GetReverse( void ) const ;
										/// 使向量反轉
			CVectorReference4FLOAT	&	Reverse( void );
										/// 檢查向量是否為零
			bool						IsZero( void ) const ;
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 + 運算子 V[4] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z, m_w + Vector.m_w]
			CVector4FLOAT				operator+( const CVectorReference4FLOAT & Vector ) const;
										/// 覆載對向量的 - 運算子 V[4] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z, m_w - Vector.m_w]
			CVector4FLOAT				operator-( const CVectorReference4FLOAT & Vector ) const;
										/// 覆載對實數的 * 運算子 V[4] = [m_x * Value, m_y * Value, m_z * Value, m_w * Value]
			CVector4FLOAT				operator*( const float Value ) const;
										/// 覆載對實數的 / 運算子 V[4] = [m_x / Value, m_y / Value, m_z / Value, m_w / Value] 
			CVector4FLOAT				operator/( const float Value ) const;

										/// 覆載 = 運算子 V = Vector
			CVectorReference4FLOAT	&	operator=( const CVectorReference4FLOAT & Vector );
										/// 覆載 = 運算子 V = Vector
			CVectorReference4FLOAT	&	operator=(const CVector4FLOAT& Vector);

										/// 覆載對向量的 += 運算子 [m_x, m_y, m_z, m_w] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z, m_w + Vector.m_w]
			CVectorReference4FLOAT	&	operator+=( const CVectorReference4FLOAT & Vector );
										/// 覆載對向量的 -= 運算子 [m_x, m_y, m_z, m_w] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z, m_w - Vector.m_w]
			CVectorReference4FLOAT	&	operator-=( const CVectorReference4FLOAT & Vector );
										/// 覆載對實數的 *= 運算子 [m_x, m_y, m_z, m_w] = [m_x * Value, m_y * Value, m_z * Value, m_w * Value]
			CVectorReference4FLOAT	&	operator*=( const float Value );
										/// 覆載對實數的 /= 運算子 [m_x, m_y, m_z, m_w] = [m_x / Value, m_y / Value, m_z / Value, m_w / Value] 
			CVectorReference4FLOAT	&	operator/=( const float Value );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 == 運算子 dx * dx + dy * dy + dz * dz + dw * dw < 0.000001
			bool						operator==( const CVectorReference4FLOAT & Vector ) const;
										/// 覆載對向量的 != 運算子 dx * dx + dy * dy + dz * dz + dw * dw >= 0.000001
			bool						operator!=( const CVectorReference4FLOAT & Vector ) const;
};


// 型別轉換用(簡化型別)
typedef  CVectorReference2FLOAT		VecR2f;
typedef  CVectorReference3FLOAT		VecR3f;
typedef  CVectorReference4FLOAT		VecR4f;

typedef  CVectorReference2FLOAT		CVectorReference2f;
typedef  CVectorReference3FLOAT		CVectorReference3f;
typedef  CVectorReference4FLOAT		CVectorReference4f;