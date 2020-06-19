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
class  CVector2FLOAT {
public:
	union {
		struct { float					m_Buffer[2];	};									///< 資料空間
		struct { float					m_x, m_y;		};									///< 空間座標
		struct { float					m_s, m_t;		};									///< 材質座標
	};
										/// 建構式
										CVector2FLOAT( float Value0 = 0.0f, float Value1 = 0.0f);
										/// 建構式 Array[2] = {m_x, m_y}
										CVector2FLOAT( float * Array );
										/// 建構式 Array[2] = {m_x, m_y}
										CVector2FLOAT(const CVectorReference2FLOAT & Vector);
										/// 拷貝建構式
										CVector2FLOAT( const CVector2FLOAT & Vector );

										/// 解構式
	virtual							   ~CVector2FLOAT();
			//////////////////////////////////////////////////////////////////////////
										/// 設定向量 
			CVector2FLOAT			&	Set(const float x, const float y );
										/// 設定向量 Array[2] = {m_x, m_y}
			CVector2FLOAT			&	Set(const float * Array );
										/// 設定向量 V = V0 + (V1-V0) * t 
			CVector2FLOAT			&	Set(const CVector2FLOAT & V0, const CVector2FLOAT & V1, const float t );
										/// 計算兩向量的 Dot (Ax * Bx + Ay * By)
			float						Dot(const CVector2FLOAT & Vector ) const;
										/// 計算兩向量的方向 逆時針 > 0 順時針 < 0
			float						Cross(const CVector2FLOAT& Vector ) const;
										/// 使兩向量交換 
			void						Swap(CVector2FLOAT& Vector );
										/// 取得兩向量的弳度
			float						GetRadian(const CVector2FLOAT& Vector) const;
										/// 取得兩向量的角度
			float						GetAngle( const CVector2FLOAT& Vector ) const;
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
			CVector2FLOAT			&	Normalize( void );
										/// 使向量反轉 
			CVector2FLOAT			&	Reverse( void );
										/// 檢查向量是否為零
			bool						IsZero( void ) const;
			////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 + 運算子 V[2] = [m_x + Vector.m_x, m_y + Vector.m_y]
			CVector2FLOAT				operator+( const CVector2FLOAT & Vector ) const;
										/// 覆載對向量的 - 運算子 V[2] = [m_x - Vector.m_x, m_y - Vector.m_y]
			CVector2FLOAT				operator-( const CVector2FLOAT & Vector ) const;
										/// 覆載對實數的 * 運算子 V[2] = [m_x * Value, m_y * Value]
			CVector2FLOAT				operator*( const float Value ) const;
										/// 覆載對實數的 / 運算子 V[2] = [m_x / Value, m_y / Value] 
			CVector2FLOAT				operator/( const float Value ) const;

										/// 覆載 = 運算子 V = Vector
			CVector2FLOAT			&	operator=( const CVector2FLOAT & Vector );
										/// 覆載 = 運算子 V = Vector
			CVector2FLOAT			&	operator=( const CVectorReference2FLOAT & Vector);

										/// 覆載對向量的 += 運算子 [m_x, m_y] = [m_x + Vector.m_x, m_y + Vector.m_y]
			CVector2FLOAT			&	operator+=( const CVector2FLOAT & Vector );
										/// 覆載對向量的 -= 運算子 [m_x, m_y] = [m_x - Vector.m_x, m_y - Vector.m_y] 
			CVector2FLOAT			&	operator-=( const CVector2FLOAT & Vector );
										/// 覆載對實數的 *= 運算子 [m_x, m_y] = [m_x * Value, m_y * Value]
			CVector2FLOAT			&	operator*=( const float Value );
										/// 覆載對實數的 /= 運算子 [m_x, m_y] = [m_x / Value, m_y / Value]
			CVector2FLOAT			&	operator/=( const float Value );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 == 運算子 dx * dx + dy * dy <  0.000001
			bool						operator==( const CVector2FLOAT & Vector ) const ;
										/// 覆載對向量的 != 運算子 dx * dx + dy * dy >= 0.000001
			bool						operator!=( const CVector2FLOAT & Vector ) const ;
};

/// 三維的向量 3x1
class CVector3FLOAT {
public:
	// 定義各分量在不同用途時使用的名稱
	union{
			struct { float               m_Buffer[3];		};				///< 資料緩衝區
			struct { float		 		 m_x, m_y, m_z;		};				///< 在空間座標時使用的名稱
			struct { float				 m_s, m_t, m_p;		};				///< 在材質座標時使用的名稱
			struct { float				 m_r, m_g, m_b;		};				///< 在顏色座標時使用的名稱
	};
										/// 建構式
										CVector3FLOAT( float Value0 = 0.0f, float Value1 = 0.0f, float Value2 = 0.0f);
										/// 建構式 Array[3] = {m_x, m_y, m_z};
										CVector3FLOAT( float * Array );
										/// 建構式 Array[3] = {m_x, m_y, m_z};
										CVector3FLOAT(const CVectorReference2FLOAT & Vector, float Value2 = 1.0f);
										/// 建構式 Array[3] = {m_x, m_y, m_z};
										CVector3FLOAT(const CVector2FLOAT & Vector, float Value2 = 1.0f);
										/// 建構式 Array[3] = {m_x, m_y, m_z};
										CVector3FLOAT( const CVectorReference3FLOAT & Vector);
										/// 拷貝建構式 										 
										CVector3FLOAT( const CVector3FLOAT & Vector );

										/// 解構式
	virtual							   ~CVector3FLOAT();

			//////////////////////////////////////////////////////////////////////////
										/// 設定向量 
			CVector3FLOAT			&	Set( const float x, const float y, const float z );
										/// 設定向量 Array[3] = {m_x, m_y, m_z};
			CVector3FLOAT			&	Set( const float * Array );
										/// 設定向量 V = V0 + (V1-V0) * t 
			CVector3FLOAT			&	Set( const CVector3FLOAT & V0, const CVector3FLOAT & V1, const float t );
										/// 計算兩向量的 Dot (Ax * Bx + Ay * By + Az * Bz)
			float						Dot( const CVector3FLOAT & Vector ) const;
										/// 計算兩向量的 Cross 
			CVector3FLOAT				Cross( const CVector3FLOAT & Vector ) const;
										/// 使兩向量交換 
			void						Swap( CVector3FLOAT & Vector );
										/// 取得兩向量的弳度
			float						GetRadian(const CVector3FLOAT & Vector) const;
										/// 取得兩向量的角度
			float						GetAngle( const CVector3FLOAT & Vector ) const;
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
			CVector3FLOAT			&	Normalize( void );
										/// 使向量反轉 
			CVector3FLOAT			&	Reverse( void );
										/// 檢查向量是否為零
			bool						IsZero( void ) const;
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 + 運算子 V[3] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z]
			CVector3FLOAT				operator+( const CVector3FLOAT & Vector ) const;
										/// 覆載對向量的 - 運算子 V[3] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z]
			CVector3FLOAT				operator-( const CVector3FLOAT & Vector ) const;
										/// 覆載對實數的 * 運算子 V[3] = [m_x * Value, m_y * Value, m_z * Value]
			CVector3FLOAT				operator*( const float Value ) const; 

										/// 覆載對實數的 / 運算子 V[3] = [m_x / Value, m_y / Value, m_z / Value] 
			CVector3FLOAT				operator/( const float Value ) const;

										/// 覆載 = 運算子 V = Vector
			CVector3FLOAT			&	operator=( const CVector3FLOAT & Vector );
										/// 覆載 = 運算子 V = Vector
			CVector3FLOAT			&	operator=( const CVectorReference3FLOAT & Vector);

										/// 覆載對向量的 += 運算子 [m_x, m_y, m_z] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z]
			CVector3FLOAT			&	operator+=( const CVector3FLOAT & Vector );
										/// 覆載對向量的 -= 運算子 [m_x, m_y, m_z] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z] 
			CVector3FLOAT			&	operator-=( const CVector3FLOAT & Vector );

										/// 覆載對實數的 *= 運算子 [m_x, m_y, m_z] = [m_x * Value, m_y * Value, m_z * Value]
			CVector3FLOAT			&	operator*=( const float Value );
										/// 覆載對實數的 /= 運算子 [m_x, m_y, m_z] = [m_x / Value, m_y / Value, m_z / Value] 
			CVector3FLOAT			&	operator/=( const float Value );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 == 運算子 dx * dx + dy * dy + dz * dz <  0.000001
			bool						operator==( const CVector3FLOAT  & Vector ) const;
										/// 覆載對向量的 != 運算子 dx * dx + dy * dy + dz * dz >= 0.000001
			bool						operator!=( const CVector3FLOAT  & Vector ) const;

};

/// 四維的向量參考 4x1
class  CVector4FLOAT {
public:
	
	// 定義各分量在不同用途時使用的名稱
	union{
			struct { float				m_Buffer[4];		};	///< 資料的緩衝區
			struct { float		 		m_x, m_y, m_z, m_w;	};	///< 在空間座標時使用的名稱
			struct { float				m_s, m_t, m_p, m_q;	};	///< 在材質座標時使用的名稱
			struct { float				m_r, m_g, m_b, m_a;	};	///< 在顏色座標時使用的名稱
	};

										/// 建構式
										CVector4FLOAT( float Value0 = 0.0f, float Value1 = 0.0f, float Value2 = 0.0f, float Value3 = 0.0f );
										/// 建構式 Array[4] = {m_x, m_y, m_z, m_w};
										CVector4FLOAT( float * Array );
										/// 建構式 Array[4] = {m_x, m_y, m_z, m_w};
										CVector4FLOAT(const CVectorReference3FLOAT & Vector, float Value3 = 1.0f);
										/// 建構式 Array[4] = {m_x, m_y, m_z, m_w}						 
										CVector4FLOAT(const CVector3FLOAT & Vector, float Value3 = 1.0f);
										/// 建構式 Array[4] = {m_x, m_y, m_z, m_w};
										CVector4FLOAT( const CVectorReference4FLOAT & Vector);
										/// 拷貝建構式 								 
										CVector4FLOAT( const CVector4FLOAT & Vector );

										/// 解構式
	virtual							   ~CVector4FLOAT();

			//////////////////////////////////////////////////////////////////////////
										/// 設定向量
			CVector4FLOAT			&	Set( const float x, const float y, const float z, const float w );
										/// 設定向量 Array[4] = {m_x, m_y, m_z, m_w}
			CVector4FLOAT			&	Set( const float * Array );									 
										/// 使兩向量交換 
			void						Swap( CVector4FLOAT  & Vector );		
			//////////////////////////////////////////////////////////////////////////
										/// 取得反轉的向量 V[4] = {-m_x, -m_y, -m_z, -m_w}
			CVector4FLOAT 				GetReverse( void ) const ;
										/// 使向量反轉
			CVector4FLOAT			&	Reverse( void );
										/// 檢查向量是否為零
			bool						IsZero( void ) const ;
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 + 運算子 V[4] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z, m_w + Vector.m_w]
			CVector4FLOAT				operator+( const CVector4FLOAT & Vector ) const;
										/// 覆載對向量的 - 運算子 V[4] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z, m_w - Vector.m_w]
			CVector4FLOAT				operator-( const CVector4FLOAT & Vector ) const;
										/// 覆載對實數的 * 運算子 V[4] = [m_x * Value, m_y * Value, m_z * Value, m_w * Value]
			CVector4FLOAT				operator*( const float Value ) const;
										/// 覆載對實數的 / 運算子 V[4] = [m_x / Value, m_y / Value, m_z / Value, m_w / Value] 
			CVector4FLOAT				operator/( const float Value ) const;

										/// 覆載 = 運算子 V = Vector
			CVector4FLOAT			&	operator=( const CVector4FLOAT & Vector );
										/// 覆載 = 運算子 V = Vector
			CVector4FLOAT			&	operator=( const CVectorReference4FLOAT& Vector);

										/// 覆載對向量的 += 運算子 [m_x, m_y, m_z, m_w] = [m_x + Vector.m_x, m_y + Vector.m_y, m_z + Vector.m_z, m_w + Vector.m_w]
			CVector4FLOAT			&	operator+=( const CVector4FLOAT & Vector );
										/// 覆載對向量的 -= 運算子 [m_x, m_y, m_z, m_w] = [m_x - Vector.m_x, m_y - Vector.m_y, m_z - Vector.m_z, m_w - Vector.m_w]
			CVector4FLOAT			&	operator-=( const CVector4FLOAT & Vector );
										/// 覆載對實數的 *= 運算子 [m_x, m_y, m_z, m_w] = [m_x * Value, m_y * Value, m_z * Value, m_w * Value]
			CVector4FLOAT			&	operator*=( const float Value );
										/// 覆載對實數的 /= 運算子 [m_x, m_y, m_z, m_w] = [m_x / Value, m_y / Value, m_z / Value, m_w / Value] 
			CVector4FLOAT			&	operator/=( const float Value );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載對向量的 == 運算子 dx * dx + dy * dy + dz * dz + dw * dw < 0.000001
			bool						operator==( const CVector4FLOAT & Vector ) const;
										/// 覆載對向量的 != 運算子 dx * dx + dy * dy + dz * dz + dw * dw >= 0.000001
			bool						operator!=( const CVector4FLOAT & Vector ) const;
};


// 型別轉換用(簡化型別)
typedef  CVector2FLOAT				Vec2f;
typedef  CVector3FLOAT				Vec3f;
typedef  CVector4FLOAT				Vec4f;

typedef  CVector2FLOAT				CVector2f;
typedef  CVector3FLOAT				CVector3f;
typedef  CVector4FLOAT				CVector4f;