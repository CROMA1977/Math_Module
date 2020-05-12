#pragma once

class CVectorReference2FLOAT;
class CVectorReference3FLOAT;
class CVectorReference4FLOAT;

class CVector2FLOAT;
class CVector3FLOAT;
class CVector4FLOAT;

class CMatrix2FLOAT;
class CMatrix3FLOAT;
class CMatrix4FLOAT;

/// 2 x 2 的矩陣
class  CMatrix2FLOAT {
public:
											 /// 對矩陣的緩衝區及各元素的定義
	union{
			struct {	float			 m_Buffer[4];	};
			struct {	float			 _M00,_M10,
										 _M01,_M11;		};
	};	
	
										/// 預定的建構式
										CMatrix2FLOAT();
										/// 以數值來指定初值的建構式(Column Matrix)
										CMatrix2FLOAT( float M00, float	M10, float	M01, float	M11 );
										/// 以陣列來指定初值的建構式 Array[4] = { _M00, _M10, _M01, _M11 }
										CMatrix2FLOAT( float * Array );
										/// 拷貝建構式
										CMatrix2FLOAT( const CMatrix2FLOAT & Matrix );
										/// 解構式
	virtual							   ~CMatrix2FLOAT();
	
			//////////////////////////////////////////////////////////////////////////
										/// 取得矩陣的緩衝區
			float					*	GetBuffer( void );
										/// 取得反矩陣
			CMatrix2FLOAT				GetInverse( void ) const;
										/// 取得轉置矩陣
			CMatrix2FLOAT				GetTranspose( void ) const;
										/// 取得 DET 的結果
			float						GetDet( void ) const;
										/// 取得 Row Vector 
			CVectorReference2FLOAT		GetRowVector( const int Index );
										/// 取得 Column Vector 
			CVectorReference2FLOAT		GetColumnVector( const int Index );
			/////////////////////////////////////////////////////////////////////////
			
										/// 以陣列指定初值 Array[4] = { _M00, _M10, _M01, _M11 }
			CMatrix2FLOAT			&	SetMatrix( const float * Array );								
										/// 設定旋轉矩陣(逆時針) M = R
			CMatrix2FLOAT			&	SetRotate( const float Angle );										
										/// 設定縮放矩陣 M = S
			CMatrix2FLOAT			&	SetScale( const float ScaleX, const float ScaleY );
										/// 載入單位矩陣
			CMatrix2FLOAT			&	LoadIdentity( void );
										/// 轉為逆矩陣
			CMatrix2FLOAT			&	Inverse( void );
										/// 轉為轉置矩陣
			CMatrix2FLOAT			&	Transpose( void );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載 + 法運算子 Mo = M + Mi
			CMatrix2FLOAT				operator+( const CMatrix2FLOAT & Matrix ) const;
										/// 覆載 - 法運算子 Mo = M - Mi
			CMatrix2FLOAT				operator-( const CMatrix2FLOAT & Matrix ) const;
										/// 覆載 * 法對實數的運算子 Mo = M * Value
			CMatrix2FLOAT				operator*( const float Value ) const; 
										/// 矩陣與向量的乘法 V = M * Vi
			CVector2FLOAT				operator*( const CVectorReference2FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector2FLOAT				operator*(const CVector2FLOAT & Vector) const;
										/// 覆載 * 法對矩陣的運算子 Mo = M * Mi
			CMatrix2FLOAT				operator*( const CMatrix2FLOAT & Matrix ) const;
										/// 覆載 / 法對實數的運算子 Mo = M / Value
			CMatrix2FLOAT				operator/( const float Value ) const;
										/// 使兩矩陣相等 M = Mi
			CMatrix2FLOAT			&	operator=( const CMatrix2FLOAT & Matrix );
										/// 覆載 += 法運算子 M = M + Mi
			CMatrix2FLOAT			&	operator+=( const CMatrix2FLOAT & Matrix );
										/// 覆載 -= 法運算子 M = M - Mi
			CMatrix2FLOAT			&	operator-=( const CMatrix2FLOAT	& Matrix );
										/// 覆載 *= 法對實數的運算子 M = M * Value
			CMatrix2FLOAT			&	operator*=( const float Value );
										/// 覆載 *= 法對矩陣的運算子 M = M * Mi
			CMatrix2FLOAT			&	operator*=( const CMatrix2FLOAT & Matrix );
										/// 覆載 /= 法對實數的運算子 M = M / Value
			CMatrix2FLOAT			&	operator/=( const float Value );
protected:
										/// 使兩矩陣相加 Mo = M + Mi
			CMatrix2FLOAT				Add( const CMatrix2FLOAT & Matrix ) const;
										/// 使兩矩陣相減 Mo = M - Mi
			CMatrix2FLOAT				Sub( const CMatrix2FLOAT & Matrix ) const;
										/// 矩陣與實數做乘法 Mo = M * Value
			CMatrix2FLOAT				Mul( const float Value ) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector2FLOAT				Mul( const CVectorReference2FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector2FLOAT				Mul(const CVector2FLOAT & Vector) const;
										/// 矩陣的乘法 Mo = M * Mi
			CMatrix2FLOAT				Mul( const CMatrix2FLOAT & Matrix ) const;
										/// 矩陣與實數做除法 M = M / Value
			CMatrix2FLOAT				Div( const float Value ) const;
										/// 使兩矩陣相等 M = Mi
			void						Assign( const CMatrix2FLOAT & Matrix );
										/// 計算兩矩陣相加並且保存數值 M = M + Mi
			void						AddAssign( const CMatrix2FLOAT & Matrix );
										/// 計算兩矩陣相減並且保存數值 M = M - Mi
			void						SubAssign( const CMatrix2FLOAT & Matrix );
										/// 計算兩矩陣相乘並且保存數值 M = M * Mi
			void						MulAssign( const CMatrix2FLOAT & Matrix );
										/// 矩陣與實數做乘法並且保存數值 M = M * Value
			void						MulAssign( const float Value );
										/// 矩陣與實數做除法並且保存數值 M = M / Value
			void						DivAssign( const float Value );
			
};

//////////////////////////////////////////////////////////////////////

/// 3 x 3 的矩陣
class  CMatrix3FLOAT {
public:	
										 /// 對矩陣的緩衝區及各元素的定義
	union{
			struct {	float			 m_Buffer[9];		};
			struct {	float			 _M00,_M10,_M20,
										 _M01,_M11,_M21,
										 _M02,_M12,_M22;	};
	};
										/// 預定的建構式
										CMatrix3FLOAT();
										/// 以數值來指定初值的建構式(Column Matrix)
										CMatrix3FLOAT(	float	M00, float	M10, float  M20,
														float	M01, float	M11, float  M21,
														float	M02, float	M12, float  M22 );
										/// 以陣列來指定初值的建構式 Array[9] = { _M00, _M10, _M20, _M01, _M11, _M21, -M02, _M12, _M22 }
										CMatrix3FLOAT( float * Array );
										/// 對 3x3 矩陣的拷貝建構式
										CMatrix3FLOAT( const CMatrix3FLOAT	& Matrix );
										/// 解構式
	virtual							   ~CMatrix3FLOAT();
			//////////////////////////////////////////////////////////////////////////
										/// 取得矩陣的緩衝區
			float					*	GetBuffer( void );
										/// 取得反矩陣
			CMatrix3FLOAT				GetInverse( void );
										/// 取得轉置矩陣
			CMatrix3FLOAT				GetTranspose( void );
										/// 取得 DET 的結果
			float						GetDet( void );
										/// 取得 Row Vector
			CVectorReference3FLOAT		GetRowVector( const int Index );
										/// 取得 Column Vector
			CVectorReference3FLOAT		GetColumnVector( const int Index );
			/////////////////////////////////////////////////////////////////////////
										/// 以陣列指定初值 Array[9] = { _M00, _M10, _M20, _M01, _M11, _M21, -M02, _M12, _M22 }
			CMatrix3FLOAT			&	SetMatrix( const float * Array );									
										/// 設定為平移矩陣 M = T2
			CMatrix3FLOAT			&	SetTranslate( const float x, const float y );
										/// 設定為平移矩陣 M = T2
			CMatrix3FLOAT			&	SetTranslate( const CVectorReference2FLOAT & Position );
										/// 設定為平移矩陣 M = T2
			CMatrix3FLOAT			&	SetTranslate(const CVector2FLOAT & Position);
										/// 設定旋轉矩陣(逆時針) M = T2 * R2 * aT2
			CMatrix3FLOAT			&	SetRotate( const float Angle, const CVectorReference2FLOAT & Center );
										/// 設定旋轉矩陣(逆時針) M = T2 * R2 * aT2
			CMatrix3FLOAT			&	SetRotate( const float Angle, const CVector2FLOAT & Center);
										/// 設定旋轉矩陣(逆時針) M = T2 * R2 * aT2
			CMatrix3FLOAT			&	SetRotate( const float Angle, const float CenterX, const float CenterY);
										/// 設定旋轉矩陣(逆時針) M = R3
			CMatrix3FLOAT			&	SetRotate( const float Angle, const CVectorReference3FLOAT & Axis );
										/// 設定旋轉矩陣(逆時針) M = R3
			CMatrix3FLOAT			&	SetRotate( const float Angle, const CVector3FLOAT & Axis);
										/// 設定旋轉矩陣(逆時針) M = R3
			CMatrix3FLOAT			&	SetRotate( const float Angle, const float AxisX, const float AxisY, const float AxisZ );
										/// 設定縮放矩陣 M = S3
			CMatrix3FLOAT			&	SetScale( const float ScaleX, const float ScaleY, const float ScaleZ = 1);
										/// 設定顏色轉換矩陣 RGB To YIQ
			CMatrix3FLOAT			&	SetColorMatrixRGB2YIQ( void );
										/// 設定顏色轉換矩陣 YIQ To RGB
			CMatrix3FLOAT			&	SetColorMatrixYIQ2RGB( void );
										/// 設定顏色轉換矩陣 RGB To YUV
			CMatrix3FLOAT			&	SetColorMatrixRGB2YUV( void );
										/// 設定顏色轉換矩陣 YUV To RGB
			CMatrix3FLOAT			&	SetColorMatrixYUV2RGB( void );
										/// 載入單位矩陣
			CMatrix3FLOAT			&	LoadIdentity( void );
										/// 計算反矩陣
			CMatrix3FLOAT			&	Inverse( void );
										/// 產生轉置矩陣
			CMatrix3FLOAT			&	Transpose( void );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載 + 法運算子 Mo = M + Mi
			CMatrix3FLOAT				operator+( const CMatrix3FLOAT & Matrix ) const;
										/// 覆載 - 法運算子 Mo = M - Mi
			CMatrix3FLOAT				operator-( const CMatrix3FLOAT & Matrix ) const;
										/// 覆載 * 法對實數的運算子 Mo = M * Value
			CMatrix3FLOAT				operator*( const float Value) const;
										/// 覆載 * 法對矩陣的運算子 Mo = M * Mi
			CMatrix3FLOAT				operator*( const CMatrix3FLOAT & Matrix ) const;
										/// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
			CVector2FLOAT				operator*( const CVectorReference2FLOAT & Vector) const;
										/// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
			CVector2FLOAT				operator*(const CVector2FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector3FLOAT				operator*( const CVectorReference3FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector3FLOAT				operator*(const CVector3FLOAT & Vector) const;
										/// 覆載 / 法對實數的運算子 Mo = M / Value
			CMatrix3FLOAT				operator/( const float Value ) const;
										/// 使兩矩陣相等 M = Mi
			CMatrix3FLOAT			&	operator=( const CMatrix3FLOAT & Matrix );
										/// 覆載 += 法運算子 M = M + Mi
			CMatrix3FLOAT			&	operator+=( const CMatrix3FLOAT & Matrix );
										/// 覆載 -= 法運算子 M = M - Mi
			CMatrix3FLOAT			&	operator-=( const CMatrix3FLOAT & Matrix );
										/// 覆載 *= 法對實數的運算子 M = M * Value
			CMatrix3FLOAT			&	operator*=( const float Value );
										/// 覆載 *= 法對矩陣的運算子 M = M * Mi
			CMatrix3FLOAT			&	operator*=( const CMatrix3FLOAT & Matrix );
										/// 覆載 /= 法對實數的運算子 M = M / Value
			CMatrix3FLOAT			&	operator/=( const float Value );
protected:
										/// 使兩矩陣相加 Mo = M + Mi
			CMatrix3FLOAT				Add( const CMatrix3FLOAT & Matrix ) const;
										/// 使兩矩陣相減 Mo = M - Mi
			CMatrix3FLOAT				Sub( const CMatrix3FLOAT & Matrix ) const;
										/// 矩陣與實數做乘法 Mo = M * Value
			CMatrix3FLOAT				Mul( const float Value ) const;
										/// 矩陣的乘法 Mo = M * Mi
			CMatrix3FLOAT				Mul( const CMatrix3FLOAT & Matrix ) const;
										/// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
			CVector2FLOAT				Mul( const CVectorReference2FLOAT & Vector) const;
										/// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
			CVector2FLOAT				Mul(const CVector2FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector3FLOAT				Mul( const CVectorReference3FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector3FLOAT				Mul(const CVector3FLOAT& Vector) const;
										/// 矩陣與實數做除法 M = M / Value
			CMatrix3FLOAT				Div( const float Value ) const;
										/// 使兩矩陣相等 M = Mi
			void						Assign( const CMatrix3FLOAT	& Matrix );
										/// 計算兩矩陣相加並且保存數值 M = M + Mi
			void						AddAssign( const CMatrix3FLOAT & Matrix );
										/// 計算兩矩陣相減並且保存數值 M = M - Mi
			void						SubAssign( const CMatrix3FLOAT & Matrix );
										/// 計算兩矩陣相乘並且保存數值 M = M * Mi
			void						MulAssign( const CMatrix3FLOAT & Matrix );
										/// 矩陣與實數做乘法並且保存數值 M = M * Value
			void						MulAssign( const float Value );
										/// 矩陣與實數做除法並且保存數值 M = M / Value
			void						DivAssign( const float Value );
};

////////////////////////////////////////////////////////////////////////////

/// 4 x 4 的矩陣
class  CMatrix4FLOAT {
public:
										/// 對矩陣的緩衝區及各元素的定義
	union{
			struct {	float			m_Buffer[16];	};
			struct {	float			_M00, _M10, _M20, _M30,
										_M01, _M11, _M21, _M31,
										_M02, _M12, _M22, _M32,
										_M03, _M13, _M23, _M33;	};
	};
										/// 預定的建構式
										CMatrix4FLOAT();
										/// 以數值來指定初值的建構式 (Column Matrix)
										CMatrix4FLOAT( float  M00, float  M10, float  M20, float  M30,
													   float  M01, float  M11, float  M21, float  M31,
													   float  M02, float  M12, float  M22, float  M32,
													   float  M03, float  M13, float  M23, float  M33);
										/// 以陣列來指定初值的建構式  Array[16] = { _M00, _M10, _M20, _M30, _M01, _M11, _M21, _M31, _M02, _M12, _M22, _M32, _M03, _M13, _M23, _M33 };
										CMatrix4FLOAT( float  * Array );
										/// 以 尤拉式 Alpha, Beta, Gamme 來初始化矩陣
										CMatrix4FLOAT( float  Alpha,	float  Beta,	float  Gamma,
													   float  OffsetX,	float  OffsetY,	float  OffsetZ );
										/// 以 四元數 Q0, Qx, Qy, Qz Position 來初始化矩陣
										CMatrix4FLOAT( float  Q0, float  Qx, float  Qy, float  Qz, const CVectorReference3FLOAT & Position );
										/// 對 4x4 矩陣的拷貝建構式
										CMatrix4FLOAT( const CMatrix4FLOAT & Matrix );		
										/// 解構式
	virtual							   ~CMatrix4FLOAT();
			//////////////////////////////////////////////////////////////////////////
										/// 取得矩陣的緩衝區
			float					*	GetBuffer( void );
										/// 取得反矩陣									
			CMatrix4FLOAT				GetInverse( void );
										/// 取得轉置矩陣
			CMatrix4FLOAT				GetTranspose( void );
										/// 計算矩陣的 Det								
			float						GetDet( void );
										/// 取得 Row Vector 
			CVectorReference4FLOAT		GetRowVector( const int Index );
										/// 取得 Column Vector 
			CVectorReference4FLOAT		GetColumnVector( const int Index );
										/// 取得轉換矩陣的軸向量
			CVectorReference3FLOAT		GetAxisVector(const int Index);
										/// 取得矩陣的姿態
			void						GetParameter( CVectorReference3FLOAT & Position, float & Rx, float & Ry, float & Rz );
			//////////////////////////////////////////////////////////////////////////
										/// 以陣列指定初值 Array[16] = { _M00, _M10, _M20, _M30, _M01, _M11, _M21, _M31, _M02, _M12, _M22, _M32, _M03, _M13, _M23, _M33 };
			CMatrix4FLOAT			&	SetMatrix( const float * Array );
										/// 以 四元數 Q0, Qx, Qy, Qz Position 初始化矩陣
			CMatrix4FLOAT			&	SetMatrix( const float Q0, const float Qx, const float Qy, const float Qz, const CVectorReference3FLOAT & Position );
										/// 以 四元數 Q0, Qx, Qy, Qz Position 初始化矩陣
			CMatrix4FLOAT			&	SetMatrix( const float Q0, const float Qx, const float Qy, const float Qz, const CVector3FLOAT& Position);
										/// 設定平移矩陣 M = T3
			CMatrix4FLOAT			&	SetTranslate( const float x, const float y, const float z );
										/// 設定平移矩陣 M = T3
			CMatrix4FLOAT			&	SetTranslate( const CVectorReference3FLOAT & Position );		
										/// 設定平移矩陣 M = T3
			CMatrix4FLOAT			&	SetTranslate(const CVector3FLOAT & Position);
										/// 設定旋轉矩陣(逆時針) M = R3
			CMatrix4FLOAT			&	SetRotate( const float Angle, const float AxisX, const float AxisY, const float AxisZ );
										/// 設定旋轉矩陣(逆時針) M = R3
			CMatrix4FLOAT			&	SetRotate( const float Angle, const CVectorReference3FLOAT & Axis );
										/// 設定旋轉矩陣(逆時針) M = R3
			CMatrix4FLOAT			&	SetRotate(const float Angle, const CVector3FLOAT & Axis);
										/// 設定旋轉矩陣(逆時針) M = T3 * R3 * aT3
			CMatrix4FLOAT			&	SetRotate( const float Angle, const CVectorReference3FLOAT & Axis, const CVectorReference3FLOAT & Center );
										/// 設定旋轉矩陣(逆時針) M = T3 * R3 * aT3
			CMatrix4FLOAT			&	SetRotate(const float Angle, const CVector3FLOAT & Axis, const CVector3FLOAT & Center);
										/// 設定縮放矩陣 M = S3
			CMatrix4FLOAT			&	SetScale( const float ScaleX, const float ScaleY, const float ScaleZ);
										/// 設定為 Ortho 投射矩陣 M = Mo
			CMatrix4FLOAT			&	Ortho( const float Left, const float Right, const float Bottom, const float Top, const float Near, const float Far );
										/// 設定為 Perspective 投射矩陣 M = Mp
			CMatrix4FLOAT			&	Perspective( const float Fovy, const float Aspect, const float zNear, const float zFar );
										/// 設定為 Frustum 投射矩陣 M = Mf
			CMatrix4FLOAT			&	Frustum( const float Left, const float Right, const float Bottom, const float Top, const float Near, const float Far );
										/// 設定為 OpenCV Camera Parameter Width, Height, Array[9]
			CMatrix4FLOAT			&   CvCameraParameter( const float Width, const float Height, const float * Array);
										/// 設定為 Lookat 攝影機矩陣
			CMatrix4FLOAT			&	LookAt( const CVectorReference3FLOAT & EyePos, const CVectorReference3FLOAT & LookPoint, const CVectorReference3FLOAT & Up );
										/// 設定為 Lookat 攝影機矩陣
			CMatrix4FLOAT			&	LookAt(const CVector3FLOAT & EyePos, const CVector3FLOAT & LookPoint, const CVector3FLOAT& Up);

										/// 載入單位矩陣
			CMatrix4FLOAT			&	LoadIdentity( void );
										/// 計算逆矩陣
			CMatrix4FLOAT			&	Inverse( void );
										/// 產生轉置矩陣
			CMatrix4FLOAT			&	Transpose( void );
			//////////////////////////////////////////////////////////////////////////
										/// 覆載 + 法運算子 Mo = M + Mi
			CMatrix4FLOAT				operator+( const CMatrix4FLOAT & Matrix ) const;
										/// 覆載 - 法運算子 Mo = M - Mi
			CMatrix4FLOAT				operator-( const CMatrix4FLOAT & Matrix ) const;
										/// 覆載 * 法對實數的運算子	Mo = M * Value
			CMatrix4FLOAT				operator*( const float Value ) const;
										/// 覆載 * 法對矩陣的運算子 Mo = M * Mi
			CMatrix4FLOAT				operator*( const CMatrix4FLOAT & Matrix ) const;
										/// 矩陣與向量的乘法 [m_x, m_y, m_z, 1] = M * [m_x, m_y, m_z, 1]
			CVector3FLOAT				operator*( const CVectorReference3FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector4FLOAT				operator*( const CVectorReference4FLOAT & Vector) const;

										/// 覆載 / 法對實數的運算子 Mo = M / Value
			CMatrix4FLOAT				operator/( const float Value ) const;
										/// 使兩矩陣相等 M = Mi
			CMatrix4FLOAT			&	operator=( const CMatrix4FLOAT & Matrix );
										/// 覆載 += 法運算子 M = M + Mi
			CMatrix4FLOAT			&	operator+=( const CMatrix4FLOAT & Matrix );
										/// 覆載 -= 法運算子 M = M - Mi
			CMatrix4FLOAT			&	operator-=( const CMatrix4FLOAT & Matrix );
										/// 覆載 *= 法運算子 M = M * Value
			CMatrix4FLOAT			&	operator*=( const float Value );
										/// 覆載 *= 法對矩陣的運算子 M = M * Mi
			CMatrix4FLOAT			&	operator*=( const CMatrix4FLOAT & Matrix );
										/// 覆載 /= 法對實數的運算子 M = M / Value
			CMatrix4FLOAT			&	operator/=( const float Value );
protected:
										/// 使兩矩陣相加 Mo = M + Mi
			CMatrix4FLOAT				Add( const CMatrix4FLOAT & Matrix ) const;
										/// 使兩矩陣相減 Mo = M - Mi
			CMatrix4FLOAT				Sub( const CMatrix4FLOAT & Matrix ) const;
										///  矩陣與實數做乘法 Mo = M * Value
			CMatrix4FLOAT				Mul( const float Value ) const;
										/// 矩陣的乘法 Mo = M * Mi
			CMatrix4FLOAT				Mul( const CMatrix4FLOAT & Matrix ) const;
										/// 矩陣與向量的乘法 [m_x, m_y, m_z, 1] = M * [m_x, m_y, m_z, 1]
			CVector3FLOAT				Mul( const CVectorReference3FLOAT & Vector) const;
										/// 矩陣與向量的乘法 V = M * Vi
			CVector4FLOAT				Mul( const CVectorReference4FLOAT & Vector) const;
										///  矩陣與實數做除法 Mo = M / Value
			CMatrix4FLOAT				Div( const float Value ) const;
										/// 使兩矩陣相等 M = Mi
			void						Assign( const CMatrix4FLOAT  & Matrix );
										/// 計算兩矩陣相加並且保存數值 M = M + Mi
			void						AddAssign( const CMatrix4FLOAT & Matrix );
										/// 計算兩矩陣相減並且保存數值 M = M - Mi
			void						SubAssign( const CMatrix4FLOAT & Matrix );
										/// 計算兩矩陣相乘並且保存數值 M = M * Mi
			void						MulAssign( const CMatrix4FLOAT & Matrix );
										///  矩陣與實數做乘法並且保存數值 M = M * Value
			void						MulAssign( const float Value );
										///  矩陣與實數做除法並且保存數值 M = M / Value
			void						DivAssign( const float Value );
};

typedef  CMatrix2FLOAT				CMatrix2f;
typedef  CMatrix3FLOAT				CMatrix3f;
typedef  CMatrix4FLOAT				CMatrix4f;