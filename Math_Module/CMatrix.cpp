#include <math.h>
#include <vector>
#include <algorithm>

#include "CMatrix.h"
#include "CVectorReference.h"
#include "CVector.h"

// 預定的建構式
CMatrix2FLOAT::CMatrix2FLOAT() {
	_M00 = 1.0f; _M01 = 0.0f;
	_M10 = 0.0f; _M11 = 1.0f;
}

// 以數值來指定初值的建構式(Column Matrix)
CMatrix2FLOAT::CMatrix2FLOAT(float M00, float M10, float M01, float M11) {
	_M00 = M00; _M01 = M01;
	_M10 = M10; _M11 = M11;
}

// 以陣列來指定初值的建構式 Array[4] = { _M00, _M10, _M01, _M11 }
CMatrix2FLOAT::CMatrix2FLOAT(float * Array) {
	if (Array == NULL) {
		_M00 = 1.0f; _M01 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f;
	}
	else {
		// _Mxx = (float) *Array 然後 Array 平移到下一個
		_M00 = *Array++;
		_M10 = *Array++;
		_M01 = *Array++;
		_M11 = *Array++;
	}
}

// 拷貝建構式
CMatrix2FLOAT::CMatrix2FLOAT(const CMatrix2FLOAT & Matrix) {
	_M00 = Matrix._M00; _M01 = Matrix._M01;
	_M10 = Matrix._M10;	_M11 = Matrix._M11;
}

// 解構式
CMatrix2FLOAT::~CMatrix2FLOAT() {
}

// 取得矩陣的緩衝區
float* CMatrix2FLOAT::GetBuffer(void) {
	return m_Buffer;
}

// 取得反矩陣
CMatrix2FLOAT CMatrix2FLOAT::GetInverse(void) const {
	float Det = GetDet();
	// 如果 Det = 0 代表沒有逆矩陣
	if (Det == 0) {
		fprintf(stderr, "\n error:Matrix2f GetInverse Det == 0!");
		return CMatrix2FLOAT(1.0f, 0.0f, 0.0f, 1.0f);
	}
	return CMatrix2FLOAT(_M11 / Det, -_M10 / Det, -_M01 / Det, _M00 / Det);
}

// 取得轉置矩陣
CMatrix2FLOAT CMatrix2FLOAT::GetTranspose(void) const {
	return CMatrix2FLOAT(_M00, _M01, _M10, _M11);
}

// 取得 DET 的結果
float CMatrix2FLOAT::GetDet(void) const {
	return _M00 * _M11 - _M01 * _M10;
}

// 取得 Row Vector 
CVectorReference2FLOAT CMatrix2FLOAT::GetRowVector(const int Index) {

	switch (Index) {
	case 0:
		return CVectorReference2FLOAT(_M00, _M01);
	case 1:
		return CVectorReference2FLOAT(_M10, _M11);
	}
	fprintf(stderr, "\n error:Matrix2f GetRowVector Over.");
	return CVectorReference2FLOAT(this->_M00, this->_M01);
}

// 取得 Column Vector 
CVectorReference2FLOAT CMatrix2FLOAT::GetColumnVector(const int Index) {
	switch (Index) {
	case 0:
		return CVectorReference2FLOAT(_M00, _M10);
	case 1:
		return CVectorReference2FLOAT(_M01, _M11);
	}
	fprintf(stderr, "\n error:Matrix2f GetColumnVector Over.");
	return CVectorReference2FLOAT(_M00, _M10);
}

// 以陣列指定初值 Array[4] = { _M00, _M10, _M01, _M11 }
CMatrix2FLOAT& CMatrix2FLOAT::SetMatrix(const float * Array) {
	if (Array == NULL) {
		fprintf(stderr, "\n error:Matrix2f SetMatrix Array == NULL!");
		_M00 = 1.0f; _M01 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f;
	}
	else {
		// _Mxx = (float) *Array 然後 Array 平移到下一個
		_M00 = *Array++;
		_M10 = *Array++;
		_M01 = *Array++;
		_M11 = *Array++;
	}
	return *this;
}

// 設定旋轉矩陣(逆時針) M = R
CMatrix2FLOAT& CMatrix2FLOAT::SetRotate(float Angle) {
	float Cos = cosf(3.1415926f / 180.0f * Angle);
	float Sin = sinf(3.1415926f / 180.0f * Angle);

	_M00 = Cos;
	_M10 = Sin;
	_M01 = -Sin;
	_M11 = Cos;

	return *this;
}

// 設定縮放矩陣 M = S
CMatrix2FLOAT& CMatrix2FLOAT::SetScale(float ScaleX, float ScaleY) {
	CMatrix2FLOAT Matrix(ScaleX, 0, 0, ScaleY);

	_M00 = ScaleX;
	_M10 = 0.0f;
	_M01 = 0.0f;
	_M11 = ScaleY;

	return *this;
}

// 載入單位矩陣
CMatrix2FLOAT& CMatrix2FLOAT::LoadIdentity(void) {
	_M00 = 1.0f; _M01 = 0.0f;
	_M10 = 0.0f; _M11 = 1.0f;
	return *this;
}

// 計算逆矩陣
CMatrix2FLOAT& CMatrix2FLOAT::Inverse(void) {
	float Det = GetDet();
	if (Det == 0) {
		fprintf(stderr, "\n error:Matrix2f Inverse Det == 0!");
		*this = CMatrix2FLOAT(1.0f, 0.0f, 0.0f, 1.0f);
	}
	else {
		*this = CMatrix2FLOAT(_M11 / Det, -_M10 / Det, -_M01 / Det, _M00 / Det);
	}
	return *this;
}

// 產生轉置矩陣
CMatrix2FLOAT& CMatrix2FLOAT::Transpose(void) {
	*this = CMatrix2FLOAT(_M00, _M01, _M10, _M11);
	return *this;
}

// 覆載 + 法運算子 Mo = M + Mi
CMatrix2FLOAT CMatrix2FLOAT::operator+(const CMatrix2FLOAT& Matrix) const {
	return Add(Matrix);
}

// 覆載 - 法運算子 Mo = M - Mi
CMatrix2FLOAT CMatrix2FLOAT::operator-(const CMatrix2FLOAT& Matrix) const {
	return Sub(Matrix);
}

// 覆載 * 法對實數的運算子 Mo = M * Value
CMatrix2FLOAT CMatrix2FLOAT::operator*(const float Value) const {
	return Mul(Value);
}

// 覆載 * 法對矩陣的運算子 Mo = M * Mi
CMatrix2FLOAT CMatrix2FLOAT::operator*(const CMatrix2FLOAT & Matrix) const {
	return Mul(Matrix);
}

// 矩陣與向量的乘法 V = M * Vi
CVector2FLOAT CMatrix2FLOAT::operator*(const CVectorReference2FLOAT & Vector) const {
	return Mul(Vector);
}

// 矩陣與向量的乘法 V = M * Vi
CVector2FLOAT CMatrix2FLOAT::operator*(const CVector2FLOAT & Vector) const {
	return Mul(Vector);
}

// 覆載 / 法對實數的運算子 Mo = M / Value
CMatrix2FLOAT CMatrix2FLOAT::operator/(const float Value) const {
	return Div(Value);
}

// 使兩矩陣相等 M = Mi
CMatrix2FLOAT& CMatrix2FLOAT::operator=(const CMatrix2FLOAT& Matrix) {
	Assign(Matrix);
	return *this;
}

// 覆載 += 法運算子 M = M + Mi
CMatrix2FLOAT& CMatrix2FLOAT::operator+=(const CMatrix2FLOAT& Matrix) {
	AddAssign(Matrix);
	return *this;
}

// 覆載 -= 法運算子 M = M - Mi
CMatrix2FLOAT& CMatrix2FLOAT::operator-=(const CMatrix2FLOAT& Matrix) {
	SubAssign(Matrix);
	return *this;
}

// 覆載 *= 法對實數的運算子 M = M * Value
CMatrix2FLOAT& CMatrix2FLOAT::operator*=(const float Value) {
	MulAssign(Value);
	return *this;
}

// 覆載 *= 法對矩陣的運算子 M = M * Mi
CMatrix2FLOAT& CMatrix2FLOAT::operator*=(const CMatrix2FLOAT & Matrix) {
	MulAssign(Matrix);
	return *this;
}

// 覆載 /= 法對實數的運算子 M = M / Value
CMatrix2FLOAT& CMatrix2FLOAT::operator/=(float Value) {
	DivAssign(Value);
	return *this;
}

// 使兩矩陣相加 Mo = M + Mi
CMatrix2FLOAT CMatrix2FLOAT::Add(const CMatrix2FLOAT & Matrix) const {
	return CMatrix2FLOAT(_M00 + Matrix._M00,
		_M10 + Matrix._M10,
		_M01 + Matrix._M01,
		_M11 + Matrix._M11);
}

// 使兩矩陣相減 Mo = M - Mi
CMatrix2FLOAT CMatrix2FLOAT::Sub(const CMatrix2FLOAT & Matrix) const {
	return CMatrix2FLOAT(_M00 - Matrix._M00,
		_M10 - Matrix._M10,
		_M01 - Matrix._M01,
		_M11 - Matrix._M11);

}

//  矩陣與實數做乘法 Mo = M * Value
CMatrix2FLOAT CMatrix2FLOAT::Mul(const float Value) const {
	return CMatrix2FLOAT(_M00 * Value,
		_M10 * Value,
		_M01 * Value,
		_M11 * Value);
}

// 矩陣的乘法 Mo = M * Mi
CMatrix2FLOAT CMatrix2FLOAT::Mul(const CMatrix2FLOAT & Matrix) const {
	return CMatrix2FLOAT(_M00 * Matrix._M00 + _M01 * Matrix._M10,
		_M10 * Matrix._M00 + _M11 * Matrix._M10,
		_M00 * Matrix._M01 + _M01 * Matrix._M11,
		_M10 * Matrix._M01 + _M11 * Matrix._M11);
}

// 矩陣與向量的乘法 V = M * Vi
CVector2FLOAT CMatrix2FLOAT::Mul(const CVectorReference2FLOAT & Vector) const {
	return CVector2FLOAT(_M00 * Vector.m_x + _M01 * Vector.m_y,
		_M10 * Vector.m_x + _M11 * Vector.m_y);
}
// 矩陣與向量的乘法 V = M * Vi
CVector2FLOAT CMatrix2FLOAT::Mul(const CVector2FLOAT & Vector) const
{
	return CVector2FLOAT(_M00 * Vector.m_x + _M01 * Vector.m_y,	_M10 * Vector.m_x + _M11 * Vector.m_y);
}

// 矩陣與實數做除法 M = M / Value
CMatrix2FLOAT CMatrix2FLOAT::Div(const float Value) const {
	if (Value == 0) {
		fprintf(stderr, "\n error:Matrix2f Div Value == 0!");
		return *this;
	}
	return CMatrix2FLOAT(_M00 / Value, _M10 / Value, _M01 / Value, _M11 / Value);
}

// 使兩矩陣相等 M = Mi
void CMatrix2FLOAT::Assign(const CMatrix2FLOAT & Matrix) {
	_M00 = Matrix._M00;
	_M10 = Matrix._M10;
	_M01 = Matrix._M01;
	_M11 = Matrix._M11;
}

// 計算兩矩陣相加並且保存數值 M = M + Mi
void CMatrix2FLOAT::AddAssign(const CMatrix2FLOAT & Matrix) {
	_M00 += Matrix._M00;
	_M10 += Matrix._M10;
	_M01 += Matrix._M01;
	_M11 += Matrix._M11;
}

// 計算兩矩陣相減並且保存數值 M = M - Mi
void CMatrix2FLOAT::SubAssign(const CMatrix2FLOAT & Matrix) {
	_M00 -= Matrix._M00;
	_M10 -= Matrix._M10;
	_M01 -= Matrix._M01;
	_M11 -= Matrix._M11;
}

// 計算兩矩陣相乘並且保存數值 M = M * Mi
void CMatrix2FLOAT::MulAssign(const CMatrix2FLOAT & Matrix) {
	*this = (*this).Mul(Matrix);
}

// 矩陣與實數做乘法並且保存數值 M = M * Value
void CMatrix2FLOAT::MulAssign(const float Value) {
	_M00 *= Value;
	_M10 *= Value;
	_M01 *= Value;
	_M11 *= Value;
}

// 矩陣與實數做除法並且保存數值 M = M / Value
void CMatrix2FLOAT::DivAssign(const float Value) {
	if (Value == 0) {
		fprintf(stderr, "\n error:Matrix2f DivAssign Value == 0!");
		return;
	}
	_M00 /= Value; _M10 /= Value; _M01 /= Value; _M11 /= Value;
}

//////////////////////////////////////////////////////////////////////



// 預定的建構式
CMatrix3FLOAT::CMatrix3FLOAT() {
	//memcpy(m_Buffer, m_Identity, sizeof(m_Identity));
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

}

// 以數值來指定初值的建構式(Column Matrix)
CMatrix3FLOAT::CMatrix3FLOAT(float M00, float M10, float M20, float M01, float M11, float M21, float M02, float M12, float M22) {
	_M00 = M00; _M10 = M10; _M20 = M20;
	_M01 = M01; _M11 = M11; _M21 = M21;
	_M02 = M02; _M12 = M12; _M22 = M22;
}

// 以陣列來指定初值的建構式 Array[9] = { _M00, _M10, _M20, _M01, _M11, _M21, -M02, _M12, _M22 }
CMatrix3FLOAT::CMatrix3FLOAT(float * Array) {
	if (Array == NULL) {
		fprintf(stderr, "\n error:Matrix3f CMatrix3FLOAT Array == NULL!");
		_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f;
		_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;
	}
	else {
		// _Mxx = (float) *Array 然後 Array 平移到下一個
		_M00 = *Array++;
		_M10 = *Array++;
		_M20 = *Array++;
		_M01 = *Array++;
		_M11 = *Array++;
		_M21 = *Array++;
		_M02 = *Array++;
		_M12 = *Array++;
		_M22 = *Array++;
	}
}

// 對 3x3 矩陣的拷貝建構式
CMatrix3FLOAT::CMatrix3FLOAT(const CMatrix3FLOAT & Matrix) {
	_M00 = Matrix._M00; _M10 = Matrix._M10; _M20 = Matrix._M20;
	_M01 = Matrix._M01; _M11 = Matrix._M11; _M21 = Matrix._M21;
	_M02 = Matrix._M02; _M12 = Matrix._M12; _M22 = Matrix._M22;
}

// 解構式
CMatrix3FLOAT::~CMatrix3FLOAT() {
}

// 取得矩陣的緩衝區
float* CMatrix3FLOAT::GetBuffer(void) {
	return m_Buffer;
}

// 取得反矩陣
CMatrix3FLOAT CMatrix3FLOAT::GetInverse(void) {
	float Det = GetDet();
	// 如果 Det = 0 代表沒有逆矩陣
	if (Det == 0) {
		fprintf(stderr, "\n error:Matrix3f GetInverse Det == 0!");
		return CMatrix3FLOAT(1.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 1.0f);
	}
	float aDet = 1.0f / Det;
	return CMatrix3FLOAT((_M11 * _M22 - _M12 * _M21) * aDet, // M00 
						-(_M10 * _M22 - _M12 * _M20) * aDet, // M10
						 (_M10 * _M21 - _M20 * _M11) * aDet, // M20
						-(_M01 * _M22 - _M21 * _M02) * aDet, // M01
						 (_M00 * _M22 - _M20 * _M02) * aDet, // M11
						-(_M00 * _M21 - _M20 * _M01) * aDet, // M21
						 (_M01 * _M12 - _M11 * _M02) * aDet, // M02
						-(_M00 * _M12 - _M10 * _M02) * aDet, // M12
						 (_M00 * _M11 - _M01 * _M10) * aDet);// M22
}

// 取得轉置矩陣
CMatrix3FLOAT CMatrix3FLOAT::GetTranspose(void) {
	return CMatrix3FLOAT(	_M00, _M01, _M02,
							_M10, _M11, _M12,
							_M20, _M21, _M22);
}

// 取得 DET 的結果
float CMatrix3FLOAT::GetDet(void) {
	return _M00 * (_M11 * _M22 - _M21 * _M12) - _M01 * (_M10 * _M22 - _M20 * _M12) + _M02 * (_M10 * _M21 - _M20 * _M11);
}

// 取得 Row Vector 
CVectorReference3FLOAT CMatrix3FLOAT::GetRowVector(const int Index) {
	switch (Index) {
	case 0:
		return CVectorReference3FLOAT(_M00, _M01, _M02);
	case 1:
		return CVectorReference3FLOAT(_M10, _M11, _M12);
	case 2:
		return CVectorReference3FLOAT(_M20, _M21, _M22);
	}
	fprintf(stderr, "\n error:Matrix3f GetRowVecotr over.");
	return CVectorReference3FLOAT(_M00, _M01, _M02);
}

// 取得 Column Vector 
CVectorReference3FLOAT CMatrix3FLOAT::GetColumnVector(const int Index) {
	switch (Index) {
	case 0:
		return CVectorReference3FLOAT(_M00, _M10, _M20);
	case 1:
		return CVectorReference3FLOAT(_M01, _M11, _M21);
	case 2:
		return CVectorReference3FLOAT(_M02, _M12, _M22);
	}
	fprintf(stderr, "\n error:Matrix3f GetColumnVector over.");
	return CVectorReference3FLOAT(_M00, _M10, _M20);
}

// 以陣列指定初值 Array[9] = { _M00, _M10, _M20, _M01, _M11, _M21, -M02, _M12, _M22 }
CMatrix3FLOAT & CMatrix3FLOAT::SetMatrix(const float * Array) {
	if (Array == NULL) {
		fprintf(stderr, "\n error:Matrix3f SetMatrix Array == NULL.");
		_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f;
		_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;
	}
	else {
		// _Mxx = (float) *Array 然後 Array 平移到下一個
		_M00 = *Array++;
		_M10 = *Array++;
		_M20 = *Array++;
		_M01 = *Array++;
		_M11 = *Array++;
		_M21 = *Array++;
		_M02 = *Array++;
		_M12 = *Array++;
		_M22 = *Array++;
	}
	return *this;
}

// 設定為平移矩陣 M = T2
CMatrix3FLOAT & CMatrix3FLOAT::SetTranslate(const float x, const float y) {
	//	    [ 1.0, 0.0,  x  ]
	// T2 = [ 0.0, 1.0,  y  ]
	//      [ 0.0, 0.0, 1.0 ]

	_M00 = 1.0f; _M01 = 0.0f; _M02 = x;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

	return *this;
}

// 設定為平移矩陣 M = T2
CMatrix3FLOAT& CMatrix3FLOAT::SetTranslate(const CVectorReference2FLOAT & Position) {
	//	    [ 1.0, 0.0,  x  ]
	// T2 = [ 0.0, 1.0,  y  ]
	//      [ 0.0, 0.0, 1.0 ]

	_M00 = 1.0f; _M01 = 0.0f; _M02 = Position.m_x;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = Position.m_y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

	return *this;
}

// 設定為平移矩陣 M = T2
CMatrix3FLOAT& CMatrix3FLOAT::SetTranslate(const CVector2FLOAT & Position) {
	//	    [ 1.0, 0.0,  x  ]
	// T2 = [ 0.0, 1.0,  y  ]
	//      [ 0.0, 0.0, 1.0 ]

	_M00 = 1.0f; _M01 = 0.0f; _M02 = Position.m_x;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = Position.m_y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

	return *this;
}
// 設定旋轉矩陣(逆時針) M = T2 * R2 * aT2
CMatrix3FLOAT& CMatrix3FLOAT::SetRotate(const float Angle, const CVectorReference2FLOAT & Center) {
	//	    [ 1.0, 0.0,  x  ]      [ Cos, -Sin, 0.0 ]       [ 1.0, 0.0,  -x ]
	// T2 = [ 0.0, 1.0,  y  ] R2 = [ Sin,  Cos, 0.0 ] aT2 = [ 0.0, 1.0,  -y ]
	//      [ 0.0, 0.0, 1.0 ]      [ 0.0,  0.0, 1.0 ]       [ 0.0, 0.0, 1.0 ]

	float x = Center.m_x;
	float y = Center.m_y;
	float Cos = cosf(3.1415926f / 180.0f * Angle);
	float Sin = sinf(3.1415926f / 180.0f * Angle);

	_M00 = Cos;  _M01 = -Sin; _M02 = Cos * -x + Sin * y + x;
	_M10 = Sin;  _M11 = Cos;  _M12 = Sin * -x + Cos * -y + y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

	return *this;
}
// 設定旋轉矩陣(逆時針) M = T2 * R2 * aT2
CMatrix3FLOAT& CMatrix3FLOAT::SetRotate(const float Angle, const CVector2FLOAT & Center) {
	//	    [ 1.0, 0.0,  x  ]      [ Cos, -Sin, 0.0 ]       [ 1.0, 0.0,  -x ]
	// T2 = [ 0.0, 1.0,  y  ] R2 = [ Sin,  Cos, 0.0 ] aT2 = [ 0.0, 1.0,  -y ]
	//      [ 0.0, 0.0, 1.0 ]      [ 0.0,  0.0, 1.0 ]       [ 0.0, 0.0, 1.0 ]

	float x = Center.m_x;
	float y = Center.m_y;
	float Cos = cosf(3.1415926f / 180.0f * Angle);
	float Sin = sinf(3.1415926f / 180.0f * Angle);

	_M00 = Cos;  _M01 = -Sin; _M02 = Cos * -x + Sin * y + x;
	_M10 = Sin;  _M11 = Cos;  _M12 = Sin * -x + Cos * -y + y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

	return *this;
}

// 設定旋轉矩陣(逆時針) M = T2 * R2 * aT2
CMatrix3FLOAT& CMatrix3FLOAT::SetRotate(const float Angle, const float CenterX, const float CenterY) {
	//	    [ 1.0, 0.0,  x  ]      [ Cos, -Sin, 0.0 ]       [ 1.0, 0.0,  -x ]
	// T2 = [ 0.0, 1.0,  y  ] R2 = [ Sin,  Cos, 0.0 ] aT2 = [ 0.0, 1.0,  -y ]
	//      [ 0.0, 0.0, 1.0 ]      [ 0.0,  0.0, 1.0 ]       [ 0.0, 0.0, 1.0 ]

	float x = CenterX;
	float y = CenterY;
	float Cos = cosf(3.1415926f / 180.0f * Angle);
	float Sin = sinf(3.1415926f / 180.0f * Angle);

	_M00 = Cos;  _M01 = -Sin; _M02 = Cos * -x + Sin * y + x;
	_M10 = Sin;  _M11 = Cos;  _M12 = Sin * -x + Cos * -y + y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;

	return *this;
}

// 設定旋轉矩陣(逆時針) M = R3
CMatrix3FLOAT& CMatrix3FLOAT::SetRotate(const float Angle, const CVectorReference3FLOAT & Axis) {
	float x, y, z;
	CVectorReference3FLOAT Vector(x, y, z);
	Vector = Axis;
	Vector.Normalize();

	float Cos = cosf(Angle * 3.1415926f / 180);
	float Sin = sinf(Angle * 3.1415926f / 180);

	_M00 = x * x * (1 - Cos) + Cos;     _M01 = x * y * (1 - Cos) - z * Sin; _M02 = x * z * (1 - Cos) + y * Sin;
	_M10 = y * x * (1 - Cos) + z * Sin; _M11 = y * y * (1 - Cos) + Cos;     _M12 = y * z * (1 - Cos) - x * Sin;
	_M20 = z * x * (1 - Cos) - y * Sin; _M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;

	return *this;
}
// 設定旋轉矩陣(逆時針) M = R3
CMatrix3FLOAT& CMatrix3FLOAT::SetRotate(const float Angle, const CVector3FLOAT & Axis)
{
	float x, y, z;
	CVectorReference3FLOAT Vector(x, y, z);
	Vector = Axis;
	Vector.Normalize();

	float Cos = cosf(Angle * 3.1415926f / 180);
	float Sin = sinf(Angle * 3.1415926f / 180);

	_M00 = x * x * (1 - Cos) + Cos;     _M01 = x * y * (1 - Cos) - z * Sin; _M02 = x * z * (1 - Cos) + y * Sin;
	_M10 = y * x * (1 - Cos) + z * Sin; _M11 = y * y * (1 - Cos) + Cos;     _M12 = y * z * (1 - Cos) - x * Sin;
	_M20 = z * x * (1 - Cos) - y * Sin; _M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;

	return *this;
}

// 設定旋轉矩陣(逆時針) M = R3
CMatrix3FLOAT& CMatrix3FLOAT::SetRotate(const float Angle, const float AxisX, const float AxisY, const float AxisZ) {
	CVector3FLOAT Axis(AxisX, AxisY, AxisZ);
	Axis.Normalize();

	float Cos = cosf(Angle * 3.1415926f / 180.0f);
	float Sin = sinf(Angle * 3.1415926f / 180.0f);
	float x = AxisX;
	float y = AxisY;
	float z = AxisZ;

	_M00 = x * x * (1 - Cos) + Cos;     _M01 = x * y * (1 - Cos) - z * Sin; _M02 = x * z * (1 - Cos) + y * Sin;
	_M10 = y * x * (1 - Cos) + z * Sin; _M11 = y * y * (1 - Cos) + Cos;     _M12 = y * z * (1 - Cos) - x * Sin;
	_M20 = z * x * (1 - Cos) - y * Sin; _M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;

	return *this;
}

// 設定縮放矩陣 M = S3
CMatrix3FLOAT& CMatrix3FLOAT::SetScale(const float ScaleX, const float ScaleY, const float ScaleZ) {
	_M00 = ScaleX; _M01 = 0.0f;   _M02 = 0.0f;
	_M10 = 0.0f;   _M11 = ScaleY; _M12 = 0.0f;
	_M20 = 0.0f;   _M21 = 0.0f;   _M22 = ScaleZ;

	return *this;
}

// 設定顏色轉換矩陣 RGB To YIQ
CMatrix3FLOAT& CMatrix3FLOAT::SetColorMatrixRGB2YIQ(void) {
	_M00 = 0.299f; _M01 = 0.587f;  _M02 = 0.114f;
	_M10 = 0.596f; _M11 = -0.274f; _M12 = -0.322f;
	_M20 = 0.211f; _M21 = -0.523f; _M22 = 0.312f;

	return *this;
}

// 設定顏色轉換矩陣 YIQ To RGB
CMatrix3FLOAT& CMatrix3FLOAT::SetColorMatrixYIQ2RGB(void) {
	_M00 = 1.000f; _M01 = 0.956f;  _M02 = 0.621f;
	_M10 = 1.000f; _M11 = -0.272f; _M12 = -0.647f;
	_M20 = 1.000f; _M21 = -1.106f; _M22 = -1.703f;

	return *this;
}

// 設定顏色轉換矩陣 RGB To YUV
CMatrix3FLOAT& CMatrix3FLOAT::SetColorMatrixRGB2YUV(void) {
	_M00 = 0.299f;  _M01 = 0.587f;  _M02 = 0.114f;
	_M10 = -0.148f; _M11 = -0.289f; _M12 = 0.437f;
	_M20 = 0.615f;  _M21 = -0.515f; _M22 = -0.100f;

	return *this;
}

// 設定顏色轉換矩陣 YUV To RGB
CMatrix3FLOAT& CMatrix3FLOAT::SetColorMatrixYUV2RGB(void) {
	_M00 = 1.000f; _M01 = 0.000f;	_M02 = 1.140f;
	_M10 = 1.000f; _M11 = -0.395f;	_M12 = -0.581f;
	_M20 = 1.000f; _M21 = 2.032f;	_M22 = 0.000f;

	return *this;
}

// 載入單位矩陣
CMatrix3FLOAT& CMatrix3FLOAT::LoadIdentity(void) {
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f;
	return *this;
}

// 計算反矩陣
CMatrix3FLOAT& CMatrix3FLOAT::Inverse(void) {
	*this = GetInverse();
	return *this;
}

// 產生轉置矩陣
CMatrix3FLOAT& CMatrix3FLOAT::Transpose(void) {
	*this = CMatrix3FLOAT(_M00, _M01, _M02,
		_M10, _M11, _M12,
		_M20, _M21, _M22);
	return *this;
}

// 覆載 + 法運算子 Mo = M + Mi
CMatrix3FLOAT CMatrix3FLOAT::operator+(const CMatrix3FLOAT & Matrix) const {
	return Add(Matrix);
}

// 覆載 - 法運算子 Mo = M - Mi
CMatrix3FLOAT CMatrix3FLOAT::operator-(const CMatrix3FLOAT & Matrix) const {
	return Sub(Matrix);
}

// 覆載 * 法對實數的運算子 Mo = M * Value
CMatrix3FLOAT CMatrix3FLOAT::operator*(const float Value) const {
	return Mul(Value);
}

// 覆載 * 法對矩陣的運算子 Mo = M * Mi
CMatrix3FLOAT CMatrix3FLOAT::operator*(const CMatrix3FLOAT & Matrix) const {
	return Mul(Matrix);
}

// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
CVector2FLOAT CMatrix3FLOAT::operator*(const CVectorReference2FLOAT & Vector) const {
	return Mul(Vector);
}

// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
CVector2FLOAT CMatrix3FLOAT::operator*(const CVector2FLOAT & Vector) const
{
	return Mul(Vector);
}

// 矩陣與向量的乘法 V = M * Vi
CVector3FLOAT CMatrix3FLOAT::operator*(const CVectorReference3FLOAT & Vector) const {
	return Mul(Vector);
}
// 矩陣與向量的乘法 V = M * Vi
CVector3FLOAT CMatrix3FLOAT::operator*(const CVector3FLOAT & Vector) const {
	return Mul(Vector);
}

// 覆載 / 法對實數的運算子 M = M / Value
CMatrix3FLOAT CMatrix3FLOAT::operator/(const float Value) const {
	return Div(Value);
}

// 使兩矩陣相等 M = Mi
CMatrix3FLOAT& CMatrix3FLOAT::operator=(const CMatrix3FLOAT & Matrix) {
	Assign(Matrix);
	return *this;
}

// 覆載 += 法運算子 M = M + Mi
CMatrix3FLOAT& CMatrix3FLOAT::operator+=(const CMatrix3FLOAT & Matrix) {
	AddAssign(Matrix);
	return *this;
}

// 覆載 -= 法運算子 M = M - Mi
CMatrix3FLOAT& CMatrix3FLOAT::operator-=(const CMatrix3FLOAT & Matrix) {
	SubAssign(Matrix);
	return *this;
}

// 覆載 *= 法對實數的運算子 M = M * Value
CMatrix3FLOAT& CMatrix3FLOAT::operator*=(const float Value) {
	MulAssign(Value);
	return *this;
}

// 覆載 *= 法對矩陣的運算子 M = M * Mi
CMatrix3FLOAT& CMatrix3FLOAT::operator*=(const CMatrix3FLOAT & Matrix) {
	MulAssign(Matrix);
	return *this;
}

// 覆載 /= 法對實數的運算子 M = M / Value
CMatrix3FLOAT& CMatrix3FLOAT::operator/=(const float Value) {
	DivAssign(Value);
	return *this;
}

// 使兩矩陣相加 Mo = M + Mi
CMatrix3FLOAT CMatrix3FLOAT::Add(const CMatrix3FLOAT & Matrix) const {
	return CMatrix3FLOAT(	_M00 + Matrix._M00, _M10 + Matrix._M10, _M20 + Matrix._M20,
							_M01 + Matrix._M01, _M11 + Matrix._M11, _M21 + Matrix._M21,
							_M02 + Matrix._M02, _M12 + Matrix._M12, _M22 + Matrix._M22);
}

// 使兩矩陣相減 Mo = M - Mi
CMatrix3FLOAT CMatrix3FLOAT::Sub(const CMatrix3FLOAT & Matrix) const {
	return CMatrix3FLOAT(	_M00 - Matrix._M00, _M10 - Matrix._M10, _M20 - Matrix._M20,
							_M01 - Matrix._M01, _M11 - Matrix._M11, _M21 - Matrix._M21,
							_M02 - Matrix._M02, _M12 - Matrix._M12, _M22 - Matrix._M22);
}

// 矩陣與實數做乘法 Mo = M * Value
CMatrix3FLOAT CMatrix3FLOAT::Mul(const float Value) const {
	return CMatrix3FLOAT(	_M00 * Value, _M10 * Value, _M20 * Value,
							_M01 * Value, _M11 * Value, _M21 * Value,
							_M02 * Value, _M12 * Value, _M22 * Value);
}

// 矩陣的乘法 Mo = M * Mi
CMatrix3FLOAT CMatrix3FLOAT::Mul(const CMatrix3FLOAT& Matrix) const {
	return CMatrix3FLOAT(	_M00 * Matrix._M00 + _M01 * Matrix._M10 + _M02 * Matrix._M20,  // M00
							_M00 * Matrix._M01 + _M01 * Matrix._M11 + _M02 * Matrix._M21,  // M01
							_M00 * Matrix._M02 + _M01 * Matrix._M12 + _M02 * Matrix._M22,  // M02
							_M10 * Matrix._M00 + _M11 * Matrix._M10 + _M12 * Matrix._M20,  // M10
							_M10 * Matrix._M01 + _M11 * Matrix._M11 + _M12 * Matrix._M21,  // M11
							_M10 * Matrix._M02 + _M11 * Matrix._M12 + _M12 * Matrix._M22,  // M12
							_M20 * Matrix._M00 + _M21 * Matrix._M10 + _M22 * Matrix._M20,  // M20
							_M20 * Matrix._M01 + _M21 * Matrix._M11 + _M22 * Matrix._M21,  // M21
							_M20 * Matrix._M02 + _M21 * Matrix._M12 + _M22 * Matrix._M22); // M22
}

// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
CVector2FLOAT CMatrix3FLOAT::Mul(const CVectorReference2FLOAT & Vector) const {
	float z = _M20 * Vector.m_x + _M21 * Vector.m_y + _M22;

	if (z == 0.0f) {
		fprintf(stderr, "\n error:Matrix3f Mul Vec2f z == 0.0f.");
		z = 1.0f;
	}

	z = 1.0f / z;
	return CVector2FLOAT( (_M00 * Vector.m_x + _M01 * Vector.m_y + _M02) * z, (_M10 * Vector.m_x + _M11 * Vector.m_y + _M12) * z);
}
// 矩陣與向量的乘法 [m_x, m_y, 1] = M * [m_x, m_y, 1]
CVector2FLOAT CMatrix3FLOAT::Mul(const CVector2FLOAT & Vector) const {
	float z = _M20 * Vector.m_x + _M21 * Vector.m_y + _M22;

	if (z == 0.0f) {
		fprintf(stderr, "\n error:Matrix3f Mul Vec2f z == 0.0f.");
		z = 1.0f;
	}

	z = 1.0f / z;
	return CVector2FLOAT((_M00 * Vector.m_x + _M01 * Vector.m_y + _M02) * z, (_M10 * Vector.m_x + _M11 * Vector.m_y + _M12) * z);
}

// 矩陣與向量的乘法 V = M * Vi
CVector3FLOAT CMatrix3FLOAT::Mul(const CVectorReference3FLOAT & Vector) const {
	return CVector3FLOAT(	_M00 * Vector.m_x + _M01 * Vector.m_y + _M02 * Vector.m_z,
							_M10 * Vector.m_x + _M11 * Vector.m_y + _M12 * Vector.m_z,
							_M20 * Vector.m_x + _M21 * Vector.m_y + _M22 * Vector.m_z);
}
// 矩陣與向量的乘法 V = M * Vi
CVector3FLOAT CMatrix3FLOAT::Mul(const CVector3FLOAT & Vector) const {
	return CVector3FLOAT(	_M00 * Vector.m_x + _M01 * Vector.m_y + _M02 * Vector.m_z,
							_M10 * Vector.m_x + _M11 * Vector.m_y + _M12 * Vector.m_z,
							_M20 * Vector.m_x + _M21 * Vector.m_y + _M22 * Vector.m_z);
}

// 矩陣與實數做除法 M = M / Value
CMatrix3FLOAT CMatrix3FLOAT::Div(const float Value) const {
	if (Value == 0) {
		fprintf(stderr, "\n error:Matrix3f Div Value == 0.0f.");
		return *this;
	}
	return CMatrix3FLOAT(	_M00 / Value, _M10 / Value, _M20 / Value,
							_M01 / Value, _M11 / Value, _M21 / Value,
							_M02 / Value, _M12 / Value, _M22 / Value);
}

// 使兩矩陣相等 M = Mi
void CMatrix3FLOAT::Assign(const CMatrix3FLOAT & Matrix) {
	_M00 = Matrix._M00; _M01 = Matrix._M01; _M02 = Matrix._M02;
	_M10 = Matrix._M10; _M11 = Matrix._M11; _M12 = Matrix._M12;
	_M20 = Matrix._M20; _M21 = Matrix._M21; _M22 = Matrix._M22;
}

// 計算兩矩陣相加並且保存數值 M = M + Mi
void CMatrix3FLOAT::AddAssign(const CMatrix3FLOAT & Matrix) {
	_M00 += Matrix._M00; _M01 += Matrix._M01; _M02 += Matrix._M02;
	_M10 += Matrix._M10; _M11 += Matrix._M11; _M12 += Matrix._M12;
	_M20 += Matrix._M20; _M21 += Matrix._M21; _M22 += Matrix._M22;
}

// 計算兩矩陣相減並且保存數值 M = M - Mi
void CMatrix3FLOAT::SubAssign(const CMatrix3FLOAT & Matrix) {
	_M00 -= Matrix._M00; _M01 -= Matrix._M01; _M02 -= Matrix._M02;
	_M10 -= Matrix._M10; _M11 -= Matrix._M11; _M12 -= Matrix._M12;
	_M20 -= Matrix._M20; _M21 -= Matrix._M21; _M22 -= Matrix._M22;
}

// 計算兩矩陣相乘並且保存數值 M = M * Mi
void CMatrix3FLOAT::MulAssign(const CMatrix3FLOAT & Matrix) {
	*this = (*this).Mul(Matrix);
}

// 矩陣與實數做乘法並且保存數值 M = M * Value
void CMatrix3FLOAT::MulAssign(const float Value) {
	_M00 *= Value; 	_M01 *= Value;	_M02 *= Value;
	_M10 *= Value; 	_M11 *= Value;	_M12 *= Value;
	_M20 *= Value;	_M21 *= Value;	_M22 *= Value;
}

void CMatrix3FLOAT::DivAssign(const float Value) {
	if (Value == 0) {
		fprintf(stderr, "\n error:Matrix3f DivAssign Value == 0.0f.");
		return;
	}
	_M00 /= Value; 	_M01 /= Value;	_M02 /= Value;
	_M10 /= Value; 	_M11 /= Value;	_M12 /= Value;
	_M20 /= Value;	_M21 /= Value;	_M22 /= Value;
}

//////////////////////////////////////////////////////////////////////////

// 預定的建構式
CMatrix4FLOAT::CMatrix4FLOAT() {
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = 0.0f;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = 0.0f;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = 0.0f;
	_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;

}

// 以數值來指定初值的建構式 (Column Matrix)
CMatrix4FLOAT::CMatrix4FLOAT(float M00, float M10, float M20, float M30, float M01, float M11, float M21, float M31, float M02, float M12, float M22, float M32, float M03, float M13, float M23, float M33) {
	_M00 = M00; _M01 = M01; _M02 = M02; _M03 = M03;
	_M10 = M10; _M11 = M11; _M12 = M12; _M13 = M13;
	_M20 = M20; _M21 = M21; _M22 = M22; _M23 = M23;
	_M30 = M30; _M31 = M31; _M32 = M32; _M33 = M33;
}

// 以陣列來指定初值的建構式  Array[16] = { _M00, _M10, _M20, _M30, _M01, _M11, _M21, _M31, _M02, _M12, _M22, _M32, _M03, _M13, _M23, _M33 };
CMatrix4FLOAT::CMatrix4FLOAT(float * Array) {
	if (Array == NULL) {
		fprintf(stderr, "\n error:Matrix4f CMatrix4FLOAT Array == NULL.");
		_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = 0.0f;
		_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = 0.0f;
		_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;
	}
	else {
		// _Mxx = (float) *Array 然後 Array 平移到下一個
		_M00 = *Array++;
		_M10 = *Array++;
		_M20 = *Array++;
		_M30 = *Array++;
		_M01 = *Array++;
		_M11 = *Array++;
		_M21 = *Array++;
		_M31 = *Array++;
		_M02 = *Array++;
		_M12 = *Array++;
		_M22 = *Array++;
		_M32 = *Array++;
		_M03 = *Array++;
		_M13 = *Array++;
		_M23 = *Array++;
		_M33 = *Array++;
	}
}

// 以 尤拉式 Alpha, Beta, Gamme 來初始化矩陣
CMatrix4FLOAT::CMatrix4FLOAT(float Alpha, float Beta, float Gamma, float OffsetX, float OffsetY, float OffsetZ) {

	float CosAlpha	= cosf(3.1415926f / 180.0f * Alpha);
	float CosBeta	= cosf(3.1415926f / 180.0f * Beta);
	float CosGamma	= cosf(3.1415926f / 180.0f * Gamma);
	float SinAlpha	= sinf(3.1415926f / 180.0f * Alpha);
	float SinBeta	= sinf(3.1415926f / 180.0f * Beta);
	float SinGamma	= sinf(3.1415926f / 180.0f * Gamma);

	_M00 = CosAlpha * CosBeta;
	_M01 = CosAlpha * SinBeta * SinGamma - SinAlpha * CosGamma;
	_M02 = CosAlpha * SinBeta * CosGamma + SinAlpha * SinGamma;
	_M03 = OffsetX;

	_M10 = SinAlpha * CosBeta;
	_M11 = SinAlpha * SinBeta * SinGamma + CosAlpha * CosGamma;
	_M12 = SinAlpha * SinBeta * CosGamma - CosAlpha * SinGamma;
	_M13 = OffsetY;

	_M20 = -SinBeta;
	_M21 = CosBeta * SinGamma;
	_M22 = CosBeta * CosGamma;
	_M23 = OffsetZ;

	_M30 = 0.0f;
	_M31 = 0.0f;
	_M32 = 0.0f;
	_M33 = 1.0f;
}

// 以 四元數 Q0, Qx, Qy, Qz Position 來初始化矩陣
CMatrix4FLOAT::CMatrix4FLOAT( const float Q0, const float Qx, const float Qy, const float Qz, const CVectorReference3FLOAT & Position) {
	_M00 = (Q0 * Q0) + (Qx * Qx) - (Qy * Qy) - (Qz * Qz);
	_M01 = 2 * ((Qx * Qy) - (Q0 * Qz));
	_M02 = 2 * ((Qx * Qz) + (Q0 * Qy));
	_M03 = Position.m_x;

	_M10 = 2 * ((Qx * Qy) + (Q0 * Qz));
	_M11 = (Q0 * Q0) - (Qx * Qx) + (Qy * Qy) - (Qz * Qz);
	_M12 = 2 * ((Qy * Qz) - (Q0 * Qx));
	_M13 = Position.m_y;

	_M20 = 2 * ((Qx * Qz) - (Q0 * Qy));
	_M21 = 2 * ((Qy * Qz) + (Q0 * Qx));
	_M22 = (Q0 * Q0) - (Qx * Qx) - (Qy * Qy) + (Qz * Qz);
	_M23 = Position.m_z;

	_M30 = 0.0f;
	_M31 = 0.0f;
	_M32 = 0.0f;
	_M33 = 1.0f;
}

// 對 4x4 矩陣的拷貝建構式
CMatrix4FLOAT::CMatrix4FLOAT(const CMatrix4FLOAT& Matrix) {
	_M00 = Matrix._M00; _M01 = Matrix._M01; _M02 = Matrix._M02; _M03 = Matrix._M03;
	_M10 = Matrix._M10; _M11 = Matrix._M11; _M12 = Matrix._M12; _M13 = Matrix._M13;
	_M20 = Matrix._M20; _M21 = Matrix._M21; _M22 = Matrix._M22; _M23 = Matrix._M23;
	_M30 = Matrix._M30; _M31 = Matrix._M31; _M32 = Matrix._M32; _M33 = Matrix._M33;
}

// 解構式
CMatrix4FLOAT::~CMatrix4FLOAT() {

}

// 取得矩陣的緩衝區
float* CMatrix4FLOAT::GetBuffer(void) {
	return m_Buffer;
}

// 取得反矩陣
CMatrix4FLOAT CMatrix4FLOAT::GetInverse(void) {
	float Det = (float)GetDet();
	// 如果 Det = 0 代表沒有逆矩陣
	if (Det == 0) {
		fprintf(stderr, "\n error:Matrix4f GetInverse Det == 0.");
		return CMatrix4FLOAT();
	}

	float aDet = 1.0f / Det;
	// 計算各區域的伴隨矩陣
	float TDetM00, TDetM01, TDetM02, TDetM03;
	float TDetM10, TDetM11, TDetM12, TDetM13;
	float TDetM20, TDetM21, TDetM22, TDetM23;
	float TDetM30, TDetM31, TDetM32, TDetM33;


	TDetM00 =	_M11 * (_M22 * _M33 - _M23 * _M32) -
				_M12 * (_M21 * _M33 - _M23 * _M31) +
				_M13 * (_M21 * _M32 - _M22 * _M31);

	TDetM01 =	_M01 * (_M22 * _M33 - _M23 * _M32) -
				_M02 * (_M21 * _M33 - _M23 * _M31) +
				_M03 * (_M21 * _M32 - _M22 * _M31);

	TDetM02 =	_M01 * (_M12 * _M33 - _M13 * _M32) -
				_M02 * (_M11 * _M33 - _M13 * _M31) +
				_M03 * (_M11 * _M32 - _M12 * _M31);

	TDetM03 =	_M01 * (_M12 * _M23 - _M13 * _M22) -
				_M02 * (_M11 * _M23 - _M13 * _M21) +
				_M03 * (_M11 * _M22 - _M12 * _M21);

	TDetM10 =	_M10 * (_M22 * _M33 - _M23 * _M32) -
				_M12 * (_M20 * _M33 - _M23 * _M30) +
				_M13 * (_M20 * _M32 - _M22 * _M30);

	TDetM11 =	_M00 * (_M22 * _M33 - _M23 * _M32) -
				_M02 * (_M20 * _M33 - _M23 * _M30) +
				_M03 * (_M20 * _M32 - _M22 * _M30);

	TDetM12 =	_M00 * (_M12 * _M33 - _M13 * _M32) -
				_M02 * (_M10 * _M33 - _M13 * _M30) +
				_M03 * (_M10 * _M32 - _M12 * _M30);

	TDetM13 =	_M00 * (_M12 * _M23 - _M13 * _M22) -
				_M02 * (_M10 * _M23 - _M13 * _M20) +
				_M03 * (_M10 * _M22 - _M12 * _M20);

	TDetM20 =	_M10 * (_M21 * _M33 - _M23 * _M31) -
				_M11 * (_M20 * _M33 - _M23 * _M30) +
				_M13 * (_M20 * _M31 - _M21 * _M30);

	TDetM21 =	_M00 * (_M21 * _M33 - _M23 * _M31) -
				_M01 * (_M20 * _M33 - _M23 * _M30) +
				_M03 * (_M20 * _M31 - _M21 * _M30);

	TDetM22 =	_M00 * (_M11 * _M33 - _M13 * _M31) -
				_M01 * (_M10 * _M33 - _M13 * _M30) +
				_M03 * (_M10 * _M31 - _M11 * _M30);

	TDetM23 =	_M00 * (_M11 * _M23 - _M13 * _M21) -
				_M01 * (_M10 * _M23 - _M13 * _M20) +
				_M03 * (_M10 * _M21 - _M11 * _M20);

	TDetM30 =	_M10 * (_M21 * _M32 - _M22 * _M31) -
				_M11 * (_M20 * _M32 - _M22 * _M30) +
				_M12 * (_M20 * _M31 - _M21 * _M30);

	TDetM31 =	_M00 * (_M21 * _M32 - _M22 * _M31) -
				_M01 * (_M20 * _M32 - _M22 * _M30) +
				_M02 * (_M20 * _M31 - _M21 * _M30);

	TDetM32 =	_M00 * (_M11 * _M32 - _M12 * _M31) -
				_M01 * (_M10 * _M32 - _M12 * _M30) +
				_M02 * (_M10 * _M31 - _M11 * _M30);

	TDetM33 =	_M00 * (_M11 * _M22 - _M12 * _M21) -
				_M01 * (_M10 * _M22 - _M12 * _M20) +
				_M02 * (_M10 * _M21 - _M11 * _M20);

	return CMatrix4FLOAT(	 TDetM00 * aDet, -TDetM10 * aDet,  TDetM20 * aDet, -TDetM30 * aDet,
							-TDetM01 * aDet,  TDetM11 * aDet, -TDetM21 * aDet,  TDetM31 * aDet,
							 TDetM02 * aDet, -TDetM12 * aDet,  TDetM22 * aDet, -TDetM32 * aDet,
							-TDetM03 * aDet,  TDetM13 * aDet, -TDetM23 * aDet,  TDetM33 * aDet);
}

// 取得轉置矩陣
CMatrix4FLOAT CMatrix4FLOAT::GetTranspose(void) {
	return CMatrix4FLOAT(	_M00, _M01, _M02, _M03,
							_M10, _M11, _M12, _M13,
							_M20, _M21, _M22, _M23,
							_M30, _M31, _M32, _M33);
}

// 計算矩陣的 Det	
float CMatrix4FLOAT::GetDet(void) {
	float DetM00, DetM01, DetM02, DetM03;

	DetM00 =	_M11 * (_M22 * _M33 - _M23 * _M32) -
				_M12 * (_M21 * _M33 - _M23 * _M31) +
				_M13 * (_M21 * _M32 - _M22 * _M31);

	DetM01 =	_M10 * (_M22 * _M33 - _M23 * _M32) -
				_M12 * (_M20 * _M33 - _M23 * _M30) +
				_M13 * (_M20 * _M32 - _M22 * _M30);

	DetM02 =	_M10 * (_M21 * _M33 - _M23 * _M31) -
				_M11 * (_M20 * _M33 - _M23 * _M30) +
				_M13 * (_M20 * _M31 - _M21 * _M30);

	DetM03 =	_M10 * (_M21 * _M32 - _M22 * _M31) -
				_M11 * (_M20 * _M32 - _M22 * _M30) +
				_M12 * (_M20 * _M31 - _M21 * _M30);

	return _M00 * DetM00 - _M01 * DetM01 + _M02 * DetM02 - _M03 * DetM03;
}

// 取得 Row Vector 
CVectorReference4FLOAT CMatrix4FLOAT::GetRowVector(int Index) {
	switch (Index) {
	case 0:
		return CVectorReference4FLOAT(_M00, _M01, _M02, _M03);
	case 1:
		return CVectorReference4FLOAT(_M10, _M11, _M12, _M13);
	case 2:
		return CVectorReference4FLOAT(_M20, _M21, _M22, _M23);
	case 3:
		return CVectorReference4FLOAT(_M30, _M31, _M32, _M33);
	}
	fprintf(stderr, "\n error:Matrix4f GetRowVector over.");
	return CVectorReference4FLOAT(_M00, _M01, _M02, _M03);
}

// 取得 Column Vector 
CVectorReference4FLOAT CMatrix4FLOAT::GetColumnVector(int Index) {
	switch (Index) {
	case 0:
		return CVectorReference4FLOAT(_M00, _M10, _M20, _M30);
	case 1:
		return CVectorReference4FLOAT(_M01, _M11, _M21, _M31);
	case 2:
		return CVectorReference4FLOAT(_M02, _M12, _M22, _M32);
	case 3:
		return CVectorReference4FLOAT(_M03, _M13, _M23, _M33);
	}
	fprintf(stderr, "\n error:Matrix4f GetColumnVector over.");
	return CVectorReference4FLOAT(_M00, _M10, _M20, _M30);
}
// 取得轉換矩陣的軸向量
CVectorReference3FLOAT CMatrix4FLOAT::GetAxisVector(const int Index)
{
	switch (Index) {
	case 0:
		return CVectorReference3FLOAT(_M00, _M10, _M20);
	case 1:
		return CVectorReference3FLOAT(_M01, _M11, _M21);
	case 2:
		return CVectorReference3FLOAT(_M02, _M12, _M22);
	case 3:
		return CVectorReference3FLOAT(_M03, _M13, _M23);
	}
	fprintf(stderr, "\n error:Matrix4f GetColumnVector over.");
	return CVectorReference3FLOAT(_M00, _M10, _M20);

}

// 取得矩陣的姿態
void CMatrix4FLOAT::GetParameter( CVectorReference3FLOAT & Position, float & Rx, float & Ry, float & Rz) {
	Position.Set(_M03, _M13, _M23);

	float GimballLock, TRx, TRy;
	Ry = -asinf(_M02);
	GimballLock = cosf(Ry);
	Ry *= (180.0f / 3.1415926f);

	if (fabs(GimballLock) > 0.005) {
		TRx = _M22 / GimballLock;
		TRy = -_M12 / GimballLock;
		Rx = atan2f(TRy, TRx) * 180.0f / 3.1415926f;
		TRx = _M00 / GimballLock;
		TRy = -_M01 / GimballLock;
		Rz = atan2f(TRy, TRx) * 180.0f / 3.1415926f;
	}
	else {
		Rx = 0;
		TRx = _M11;
		TRy = _M10;
		Rz = atan2f(TRy, TRx) * (180.0f * 3.1415926f);
	}
}

// 以陣列指定初值 Array[16] = { _M00, _M10, _M20, _M30, _M01, _M11, _M21, _M31, _M02, _M12, _M22, _M32, _M03, _M13, _M23, _M33 };
CMatrix4FLOAT & CMatrix4FLOAT::SetMatrix(const float * Array) {
	if (Array == NULL) {
		fprintf(stderr, "\n error:Matrix4f SetMatrix Array == NULL.");
		_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = 0.0f;
		_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = 0.0f;
		_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;
	}
	else {
		// _Mxx = (float) *Array 然後 Array 平移到下一個
		_M00 = *Array++;
		_M10 = *Array++;
		_M20 = *Array++;
		_M30 = *Array++;
		_M01 = *Array++;
		_M11 = *Array++;
		_M21 = *Array++;
		_M31 = *Array++;
		_M02 = *Array++;
		_M12 = *Array++;
		_M22 = *Array++;
		_M32 = *Array++;
		_M03 = *Array++;
		_M13 = *Array++;
		_M23 = *Array++;
		_M33 = *Array++;
	}
	return *this;
}

// 以 四元數 Q0, Qx, Qy, Qz Position 初始化矩陣
CMatrix4FLOAT & CMatrix4FLOAT::SetMatrix(const float Q0, const float Qx, const float Qy, const float Qz, const CVectorReference3FLOAT & Position) {
	_M00 = (Q0 * Q0) + (Qx * Qx) - (Qy * Qy) - (Qz * Qz);
	_M01 = 2 * ((Qx * Qy) - (Q0 * Qz));
	_M02 = 2 * ((Qx * Qz) + (Q0 * Qy));
	_M03 = Position.m_x;

	_M10 = 2 * ((Qx * Qy) + (Q0 * Qz));
	_M11 = (Q0 * Q0) - (Qx * Qx) + (Qy * Qy) - (Qz * Qz);
	_M12 = 2 * ((Qy * Qz) - (Q0 * Qx));
	_M13 = Position.m_y;

	_M20 = 2 * ((Qx * Qz) - (Q0 * Qy));
	_M21 = 2 * ((Qy * Qz) + (Q0 * Qx));
	_M22 = (Q0 * Q0) - (Qx * Qx) - (Qy * Qy) + (Qz * Qz);
	_M23 = Position.m_z;

	_M30 = 0.0f;
	_M31 = 0.0f;
	_M32 = 0.0f;
	_M33 = 1.0f;
	return *this;
}
// 以 四元數 Q0, Qx, Qy, Qz Position 初始化矩陣
CMatrix4FLOAT& CMatrix4FLOAT::SetMatrix(const float Q0, const float Qx, const float Qy, const float Qz, const CVector3FLOAT& Position)
{
	_M00 = (Q0 * Q0) + (Qx * Qx) - (Qy * Qy) - (Qz * Qz);
	_M01 = 2 * ((Qx * Qy) - (Q0 * Qz));
	_M02 = 2 * ((Qx * Qz) + (Q0 * Qy));
	_M03 = Position.m_x;

	_M10 = 2 * ((Qx * Qy) + (Q0 * Qz));
	_M11 = (Q0 * Q0) - (Qx * Qx) + (Qy * Qy) - (Qz * Qz);
	_M12 = 2 * ((Qy * Qz) - (Q0 * Qx));
	_M13 = Position.m_y;

	_M20 = 2 * ((Qx * Qz) - (Q0 * Qy));
	_M21 = 2 * ((Qy * Qz) + (Q0 * Qx));
	_M22 = (Q0 * Q0) - (Qx * Qx) - (Qy * Qy) + (Qz * Qz);
	_M23 = Position.m_z;

	_M30 = 0.0f;
	_M31 = 0.0f;
	_M32 = 0.0f;
	_M33 = 1.0f;
	return *this;
}

// 設定平移矩陣 M = T3
CMatrix4FLOAT & CMatrix4FLOAT::SetTranslate(const float x, const float y, const float z) {
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = x;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = z;
	_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;

	return *this;
}

// 設定平移矩陣 M = T3
CMatrix4FLOAT & CMatrix4FLOAT::SetTranslate(const CVectorReference3FLOAT & Position) {
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = Position.m_x;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = Position.m_y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = Position.m_z;
	_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;

	return *this;
}
// 設定平移矩陣 M = T3
CMatrix4FLOAT& CMatrix4FLOAT::SetTranslate(const CVector3FLOAT & Position) {
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = Position.m_x;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = Position.m_y;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = Position.m_z;
	_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;

	return *this;
}

// 設定旋轉矩陣(逆時針) M = R3
CMatrix4FLOAT & CMatrix4FLOAT::SetRotate(const float Angle, const float AxisX, const float AxisY, const float AxisZ) {
	float x, y, z;
	CVectorReference3FLOAT Axis(x, y, z);
	Axis.Set(AxisX, AxisY, AxisZ);
	Axis.Normalize();

	float Cos = cosf(Angle * 3.1415926f / 180.0f);
	float Sin = sinf(Angle * 3.1415926f / 180.0f);

	_M00 = x * x * (1 - Cos) + Cos;		_M01 = x * y * (1 - Cos) - z * Sin; _M02 = x * z * (1 - Cos) + y * Sin; _M03 = 0.0f;
	_M10 = y * x * (1 - Cos) + z * Sin; _M11 = y * y * (1 - Cos) + Cos;		_M12 = y * z * (1 - Cos) - x * Sin; _M13 = 0.0f;
	_M20 = z * x * (1 - Cos) - y * Sin; _M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;		_M23 = 0.0f;
	_M30 = 0.0f;						_M31 = 0.0f;						_M32 = 0.0f;						_M33 = 1.0f;

	return *this;
}

// 設定旋轉矩陣(逆時針) M = R3
CMatrix4FLOAT & CMatrix4FLOAT::SetRotate(const float Angle, const CVectorReference3FLOAT & Axis) {
	float x, y, z;
	CVectorReference3FLOAT Vector(x, y, z);
	Vector = Axis;
	Vector.Normalize();

	float Cos = (float)cos(Angle * 3.1415926f / 180.0f);
	float Sin = (float)sin(Angle * 3.1415926f / 180.0f);

	_M00 = x * x * (1 - Cos) + Cos;		_M01 = x * y * (1 - Cos) - z * Sin; _M02 = x * z * (1 - Cos) + y * Sin; _M03 = 0.0f;
	_M10 = y * x * (1 - Cos) + z * Sin; _M11 = y * y * (1 - Cos) + Cos;		_M12 = y * z * (1 - Cos) - x * Sin; _M13 = 0.0f;
	_M20 = z * x * (1 - Cos) - y * Sin; _M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;		_M23 = 0.0f;
	_M30 = 0.0f;						_M31 = 0.0f;						_M32 = 0.0f;						_M33 = 1.0f;

	return *this;
}
// 設定旋轉矩陣(逆時針) M = R3
CMatrix4FLOAT& CMatrix4FLOAT::SetRotate(const float Angle, const CVector3FLOAT & Axis) {
	float x, y, z;
	CVectorReference3FLOAT Vector(x, y, z);
	Vector = Axis;
	Vector.Normalize();

	float Cos = (float)cos(Angle * 3.1415926f / 180.0f);
	float Sin = (float)sin(Angle * 3.1415926f / 180.0f);

	_M00 = x * x * (1 - Cos) + Cos;		_M01 = x * y * (1 - Cos) - z * Sin; _M02 = x * z * (1 - Cos) + y * Sin; _M03 = 0.0f;
	_M10 = y * x * (1 - Cos) + z * Sin; _M11 = y * y * (1 - Cos) + Cos;		_M12 = y * z * (1 - Cos) - x * Sin; _M13 = 0.0f;
	_M20 = z * x * (1 - Cos) - y * Sin; _M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;		_M23 = 0.0f;
	_M30 = 0.0f;						_M31 = 0.0f;						_M32 = 0.0f;						_M33 = 1.0f;

	return *this;
}

// 設定旋轉矩陣(逆時針) M = T3 * R3 * aT3
CMatrix4FLOAT & CMatrix4FLOAT::SetRotate(const float Angle, const CVectorReference3FLOAT & Axis, const CVectorReference3FLOAT & Center) {
	//		[ 1.0, 0.0, 0.0, Center.m_x ]		[ 1.0, 0.0, 0.0, -Center.m_x ]
	//		[ 0.0, 1.0, 0.0, Center.m_y ]		[ 0.0, 1.0, 0.0, -Center.m_y ]
	// T3 = [ 0.0, 0.0, 1.0, Center.m_z ] aT3 = [ 0.0, 0.0, 1.0, -Center.m_z ]
	//		[ 0.0, 0.0, 0.0, 1.0		]		[ 0.0, 0.0, 0.0,  1.0		 ]
	//		[ x * x * (1 - Cos) + Cos,			x * y * (1 - Cos) - z * Sin,	x * z * (1 - Cos) + y * Sin,	0.0 ]
	//		[ y * x * (1 - Cos) + z * Sin,		y * y * (1 - Cos) + Cos,		y * z * (1 - Cos) - x * Sin,	0.0 ]
	// R3 = [ z * x * (1 - Cos) - y * Sin,		z * y * (1 - Cos) + x * Sin,	z * z * (1 - Cos) + Cos,		0.0 ]
	// 		[ 0.0,								0.0,							0.0,							1.0 ]
	//  M = T3 * R3 * aT3;

	float x, y, z;
	CVectorReference3FLOAT Vector(x, y, z);
	Vector = Axis;
	Vector.Normalize();

	float Cos = cosf(Angle * 3.1415926f / 180.0f);
	float Sin = sinf(Angle * 3.1415926f / 180.0f);
	float CosX = Cos * x; float CosY = Cos * y; float CosZ = Cos * z;
	float SinX = Sin * x; float SinY = Sin * y; float SinZ = Sin * z;
	float one_minus_Cos = 1.0f - Cos;
	float xx = x * x * one_minus_Cos; float xy = x * y * one_minus_Cos; float xz = x * z * one_minus_Cos;
	float yy = y * y * one_minus_Cos; float yz = y * z * one_minus_Cos;
	float zz = z * z * one_minus_Cos;

	_M00 = xx + Cos;	_M01 = xy - SinZ;	_M02 = xz + SinY;	_M03 = _M00 * -Center.m_x + _M01 * -Center.m_y + _M02 * -Center.m_z + Center.m_x;
	_M10 = xy + SinZ;	_M11 = yy + Cos;	_M12 = yz - SinX;	_M13 = _M10 * -Center.m_x + _M11 * -Center.m_y + _M12 * -Center.m_z + Center.m_y;
	_M20 = xz - SinY;	_M21 = yz + SinX;	_M22 = zz + Cos;	_M23 = _M20 * -Center.m_x + _M21 * -Center.m_y + _M22 * -Center.m_z + Center.m_z;
	_M30 = 0.0f;		_M31 = 0.0f;		_M32 = 0.0f;		_M33 = 1.0f;

	//_M00 = x * x * (1 - Cos) + Cos;		_M01 = x * y * (1 - Cos) - z * Sin;	_M02 = x * z * (1 - Cos) + y * Sin;	_M03 = (x * x * (1 - Cos) + Cos) * -Center.m_x + (x * y * (1 - Cos) - z * Sin) * -Center.m_y + (x * z * (1 - Cos) + y * Sin ) * -Center.m_z + Center.m_x;
	//_M10 = y * x * (1 - Cos) + z * Sin;	_M11 = y * y * (1 - Cos) + Cos;		_M12 = y * z * (1 - Cos) - x * Sin; _M13 = (y * x * (1 - Cos) + z * Sin) * -Center.m_x + (y * y * (1 - Cos) + Cos) * -Center.m_y + (y * z * (1 - Cos) - x * Sin) * -Center.m_z  + Center.m_y;
	//_M20 = z * x * (1 - Cos) - y * Sin;	_M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;		_M23 = (z * x * (1 - Cos) - y * Sin) * -Center.m_x + (z * y * (1 - Cos) + x * Sin) * -Center.m_y + (z * z * (1 - Cos) + Cos) * -Center.m_z  + Center.m_z;
	//_M30 = 0.0f;							_M31 = 0.0f;						_M32 = 0.0f;						_M33 = 1.0f;
	return *this;
}
// 設定旋轉矩陣(逆時針) M = T3 * R3 * aT3
CMatrix4FLOAT& CMatrix4FLOAT::SetRotate(const float Angle, const CVector3FLOAT & Axis, const CVector3FLOAT & Center) {
	//		[ 1.0, 0.0, 0.0, Center.m_x ]		[ 1.0, 0.0, 0.0, -Center.m_x ]
	//		[ 0.0, 1.0, 0.0, Center.m_y ]		[ 0.0, 1.0, 0.0, -Center.m_y ]
	// T3 = [ 0.0, 0.0, 1.0, Center.m_z ] aT3 = [ 0.0, 0.0, 1.0, -Center.m_z ]
	//		[ 0.0, 0.0, 0.0, 1.0		]		[ 0.0, 0.0, 0.0,  1.0		 ]
	//		[ x * x * (1 - Cos) + Cos,			x * y * (1 - Cos) - z * Sin,	x * z * (1 - Cos) + y * Sin,	0.0 ]
	//		[ y * x * (1 - Cos) + z * Sin,		y * y * (1 - Cos) + Cos,		y * z * (1 - Cos) - x * Sin,	0.0 ]
	// R3 = [ z * x * (1 - Cos) - y * Sin,		z * y * (1 - Cos) + x * Sin,	z * z * (1 - Cos) + Cos,		0.0 ]
	// 		[ 0.0,								0.0,							0.0,							1.0 ]
	//  M = T3 * R3 * aT3;

	float x, y, z;
	CVectorReference3FLOAT Vector(x, y, z);
	Vector = Axis;
	Vector.Normalize();

	float Cos = cosf(Angle * 3.1415926f / 180.0f);
	float Sin = sinf(Angle * 3.1415926f / 180.0f);
	float CosX = Cos * x; float CosY = Cos * y; float CosZ = Cos * z;
	float SinX = Sin * x; float SinY = Sin * y; float SinZ = Sin * z;
	float one_minus_Cos = 1.0f - Cos;
	float xx = x * x * one_minus_Cos; float xy = x * y * one_minus_Cos; float xz = x * z * one_minus_Cos;
	float yy = y * y * one_minus_Cos; float yz = y * z * one_minus_Cos;
	float zz = z * z * one_minus_Cos;

	_M00 = xx + Cos;	_M01 = xy - SinZ;	_M02 = xz + SinY;	_M03 = _M00 * -Center.m_x + _M01 * -Center.m_y + _M02 * -Center.m_z + Center.m_x;
	_M10 = xy + SinZ;	_M11 = yy + Cos;	_M12 = yz - SinX;	_M13 = _M10 * -Center.m_x + _M11 * -Center.m_y + _M12 * -Center.m_z + Center.m_y;
	_M20 = xz - SinY;	_M21 = yz + SinX;	_M22 = zz + Cos;	_M23 = _M20 * -Center.m_x + _M21 * -Center.m_y + _M22 * -Center.m_z + Center.m_z;
	_M30 = 0.0f;		_M31 = 0.0f;		_M32 = 0.0f;		_M33 = 1.0f;

	//_M00 = x * x * (1 - Cos) + Cos;		_M01 = x * y * (1 - Cos) - z * Sin;	_M02 = x * z * (1 - Cos) + y * Sin;	_M03 = (x * x * (1 - Cos) + Cos) * -Center.m_x + (x * y * (1 - Cos) - z * Sin) * -Center.m_y + (x * z * (1 - Cos) + y * Sin ) * -Center.m_z + Center.m_x;
	//_M10 = y * x * (1 - Cos) + z * Sin;	_M11 = y * y * (1 - Cos) + Cos;		_M12 = y * z * (1 - Cos) - x * Sin; _M13 = (y * x * (1 - Cos) + z * Sin) * -Center.m_x + (y * y * (1 - Cos) + Cos) * -Center.m_y + (y * z * (1 - Cos) - x * Sin) * -Center.m_z  + Center.m_y;
	//_M20 = z * x * (1 - Cos) - y * Sin;	_M21 = z * y * (1 - Cos) + x * Sin; _M22 = z * z * (1 - Cos) + Cos;		_M23 = (z * x * (1 - Cos) - y * Sin) * -Center.m_x + (z * y * (1 - Cos) + x * Sin) * -Center.m_y + (z * z * (1 - Cos) + Cos) * -Center.m_z  + Center.m_z;
	//_M30 = 0.0f;							_M31 = 0.0f;						_M32 = 0.0f;						_M33 = 1.0f;
	return *this;
}
// 設定縮放矩陣 M = S3
CMatrix4FLOAT& CMatrix4FLOAT::SetScale(const float ScaleX, const float ScaleY, const float ScaleZ) {
	_M00 = ScaleX;	_M01 = 0.0f;	_M02 = 0.0f;	_M03 = 0.0f;
	_M10 = 0.0f;	_M11 = ScaleY;	_M12 = 0.0f;	_M13 = 0.0f;
	_M20 = 0.0f;	_M21 = 0.0f;	_M22 = ScaleZ;	_M23 = 0.0f;
	_M30 = 0.0f;	_M31 = 0.0f;	_M32 = 0.0f;	_M33 = 1.0f;

	return *this;
}

// 設定為 Ortho 投射矩陣 M = Mo
CMatrix4FLOAT & CMatrix4FLOAT::Ortho(float Left, float Right, float Bottom, float Top, float Near, float Far) {
	_M00 = (2.0f / (Right - Left));	_M01 = 0.0f; 					_M02 = 0.0f; 					_M03 = -(Right + Left) / (Right - Left);
	_M10 = 0.0f; 					_M11 = (2.0f / (Top - Bottom));	_M12 = 0.0f; 					_M13 = -(Top + Bottom) / (Top - Bottom);
	_M20 = 0.0f; 					_M21 = 0.0f; 					_M22 = -2.0f / (Far - Near);	_M23 = (Far + Near) / (Far - Near);
	_M30 = 0.0f; 					_M31 = 0.0f; 					_M32 = 0.0f; 					_M33 = 1.0f;

	return *this;
}

// 設定為 Perspective 投射矩陣 M = Mp
CMatrix4FLOAT & CMatrix4FLOAT::Perspective(const float Fovy, const float Aspect, const float zNear, const float zFar) {
	float f = 1.0f / tanf(Fovy / 2.0f * 3.1415926f / 180.0f);
	_M00 = f / Aspect;	_M01 = 0.0f;	_M02 = 0.0f;								_M03 = 0.0f;
	_M10 = 0.0f;		_M11 = f;		_M12 = 0.0f;								_M13 = 0.0f;
	_M20 = 0.0f;		_M21 = 0.0f;	_M22 = (zNear + zFar) / (zNear - zFar);		_M23 = 2.0f * zFar * zNear / (zNear - zFar);
	_M30 = 0.0f;		_M31 = 0.0f;	_M32 = -1.0f;								_M33 = 0.0f;

	return *this;
}

// 設定為 Frustum 投射矩陣 M = Mf
CMatrix4FLOAT & CMatrix4FLOAT::Frustum(float Left, float Right, float Bottom, float Top, float Near, float Far) {
	_M00 = 2.0f * Near / (Right - Left);	_M01 = 0.0f; 							_M02 = (Right + Left) / (Right - Left);	_M03 = 0.0f;
	_M10 = 0.0f; 						_M11 = 2.0f * Near / (Top - Bottom);	_M12 = (Top + Bottom) / (Top - Bottom);	_M13 = 0.0f;
	_M20 = 0.0f; 						_M21 = 0.0f; 							_M22 = -(Far + Near) / (Far - Near);	_M23 = -2.0f * Far * Near / (Far - Near);
	_M30 = 0.0f; 						_M31 = 0.0f; 							_M32 = -1.0f;									_M33 = 0.0f;

	return *this;
}

// 設定為 OpenCV Camera Parameter Width, Height, Array[9]
CMatrix4FLOAT & CMatrix4FLOAT::CvCameraParameter(const float Width, const float Height, const float * Array) {
	float Near = 1.0f;
	float Far = Array[0] * 10.0f;

	_M00 = 2.0f / Width * -Array[0];	_M01 = 2.0f / Width * Array[3];	_M02 = 2.0f / Width * Array[6];			_M03 = 0.0f;
	_M10 = 2.0f / Height * Array[1];	_M11 = 2.0f / Height * -Array[4];	_M12 = 2.0f / Height * Array[7];			_M13 = 0.0f;
	_M20 = Width / 2.0f - Array[2];		_M21 = Array[5] - Height / 2.0f;	_M22 = (Near + Far) / (Near - Far);			_M23 = -1.0f;
	_M30 = 0.0f;						_M31 = 0.0f;						_M32 = (2.0f * Near * Far) / (Near - Far);	_M33 = 0.0f;

	return *this;
}

// 設定為 Lookat 攝影機矩陣
CMatrix4FLOAT& CMatrix4FLOAT::LookAt(const CVectorReference3FLOAT & EyePos, const CVectorReference3FLOAT & LookPoint, const CVectorReference3FLOAT & Up) {
	// 以 EyePos 為中心求正交的三軸 F = -Z 軸 U = Y 軸 R = X 軸
	CVector3FLOAT F = (LookPoint - EyePos).GetNormalize();
	CVector3FLOAT U = Up.GetNormalize();
	CVector3FLOAT R = F.Cross(U);
	U = R.Cross(F);

	//		[  R.m_x,  R.m_y,  R.m_z, 0.0 ]			[ 1.0, 0.0, 0.0, -EyePos.m_x ]
	//		[  U.m_x,  U.m_y,  U.m_z, 0.0 ]			[ 0.0, 1.0, 0.0, -EyePos.m_y ]
	// Ml =	[ -F.m_x, -F.m_y, -F.m_z, 0.0 ]  Mt =	[ 0.0, 0.0, 1.0, -EyePos.m_z ]
	//		[  0.0,    0.0,    0.0,   1.0 ]			[ 0.0, 0.0, 0.0,       1.0   ]
	//
	// LookAt = Ml * Mt;

	_M00 = R.m_x;	_M01 = R.m_y;	_M02 = R.m_z;	_M03 = R.m_x * -EyePos.m_x + R.m_y * -EyePos.m_y + R.m_z * -EyePos.m_z;
	_M10 = U.m_x;	_M11 = U.m_y;	_M12 = U.m_z;	_M13 = U.m_x * -EyePos.m_x + U.m_y * -EyePos.m_y + U.m_z * -EyePos.m_z;
	_M20 = -F.m_x;	_M21 = -F.m_y;	_M22 = -F.m_z;	_M23 = -F.m_x * -EyePos.m_x + -F.m_y * -EyePos.m_y + -F.m_z * -EyePos.m_z;
	_M30 = 0.0f;	_M31 = 0.0f;	_M32 = 0.0f;	_M33 = 1.0f;

	return *this;
}
// 設定為 Lookat 攝影機矩陣
CMatrix4FLOAT & CMatrix4FLOAT::LookAt(const CVector3FLOAT & EyePos, const CVector3FLOAT & LookPoint, const CVector3FLOAT& Up) {
	// 以 EyePos 為中心求正交的三軸 F = -Z 軸 U = Y 軸 R = X 軸
	CVector3FLOAT F = (LookPoint - EyePos).GetNormalize();
	CVector3FLOAT U = Up.GetNormalize();
	CVector3FLOAT R = F.Cross(U);
	U = R.Cross(F);

	//		[  R.m_x,  R.m_y,  R.m_z, 0.0 ]			[ 1.0, 0.0, 0.0, -EyePos.m_x ]
	//		[  U.m_x,  U.m_y,  U.m_z, 0.0 ]			[ 0.0, 1.0, 0.0, -EyePos.m_y ]
	// Ml =	[ -F.m_x, -F.m_y, -F.m_z, 0.0 ]  Mt =	[ 0.0, 0.0, 1.0, -EyePos.m_z ]
	//		[  0.0,    0.0,    0.0,   1.0 ]			[ 0.0, 0.0, 0.0,       1.0   ]
	//
	// LookAt = Ml * Mt;

	_M00 =  R.m_x;	_M01 =  R.m_y;	_M02 =  R.m_z;	_M03 =  R.m_x * -EyePos.m_x +  R.m_y * -EyePos.m_y +  R.m_z * -EyePos.m_z;
	_M10 =  U.m_x;	_M11 =  U.m_y;	_M12 =  U.m_z;	_M13 =  U.m_x * -EyePos.m_x +  U.m_y * -EyePos.m_y +  U.m_z * -EyePos.m_z;
	_M20 = -F.m_x;	_M21 = -F.m_y;	_M22 = -F.m_z;	_M23 = -F.m_x * -EyePos.m_x + -F.m_y * -EyePos.m_y + -F.m_z * -EyePos.m_z;
	_M30 = 0.0f;	_M31 = 0.0f;	_M32 = 0.0f;	_M33 = 1.0f;

	return *this;
}

// 載入單位矩陣
CMatrix4FLOAT& CMatrix4FLOAT::LoadIdentity(void) {
	_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = 0.0f;
	_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = 0.0f;
	_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = 0.0f;
	_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;
	return *this;
}

// 計算逆矩陣
CMatrix4FLOAT& CMatrix4FLOAT::Inverse(void) {
	float Det = (float)GetDet();
	// 如果 Det = 0 代表沒有逆矩陣
	if (Det == 0) {
		fprintf(stderr, "\n error:Matrix4f Inverse Det == 0!.");
		_M00 = 1.0f; _M01 = 0.0f; _M02 = 0.0f; _M03 = 0.0f;
		_M10 = 0.0f; _M11 = 1.0f; _M12 = 0.0f; _M13 = 0.0f;
		_M20 = 0.0f; _M21 = 0.0f; _M22 = 1.0f; _M23 = 0.0f;
		_M30 = 0.0f; _M31 = 0.0f; _M32 = 0.0f; _M33 = 1.0f;
		return *this;
	}

	float aDet = 1.0f / Det;
	// 計算各區域的伴隨矩陣
	float TDetM00, TDetM01, TDetM02, TDetM03;
	float TDetM10, TDetM11, TDetM12, TDetM13;
	float TDetM20, TDetM21, TDetM22, TDetM23;
	float TDetM30, TDetM31, TDetM32, TDetM33;

	TDetM00 = _M11 * (_M22 * _M33 - _M23 * _M32) -
		_M12 * (_M21 * _M33 - _M23 * _M31) +
		_M13 * (_M21 * _M32 - _M22 * _M31);

	TDetM01 = _M01 * (_M22 * _M33 - _M23 * _M32) -
		_M02 * (_M21 * _M33 - _M23 * _M31) +
		_M03 * (_M21 * _M32 - _M22 * _M31);

	TDetM02 = _M01 * (_M12 * _M33 - _M13 * _M32) -
		_M02 * (_M11 * _M33 - _M13 * _M31) +
		_M03 * (_M11 * _M32 - _M12 * _M31);

	TDetM03 = _M01 * (_M12 * _M23 - _M13 * _M22) -
		_M02 * (_M11 * _M23 - _M13 * _M21) +
		_M03 * (_M11 * _M22 - _M12 * _M21);

	TDetM10 = _M10 * (_M22 * _M33 - _M23 * _M32) -
		_M12 * (_M20 * _M33 - _M23 * _M30) +
		_M13 * (_M20 * _M32 - _M22 * _M30);

	TDetM11 = _M00 * (_M22 * _M33 - _M23 * _M32) -
		_M02 * (_M20 * _M33 - _M23 * _M30) +
		_M03 * (_M20 * _M32 - _M22 * _M30);

	TDetM12 = _M00 * (_M12 * _M33 - _M13 * _M32) -
		_M02 * (_M10 * _M33 - _M13 * _M30) +
		_M03 * (_M10 * _M32 - _M12 * _M30);

	TDetM13 = _M00 * (_M12 * _M23 - _M13 * _M22) -
		_M02 * (_M10 * _M23 - _M13 * _M20) +
		_M03 * (_M10 * _M22 - _M12 * _M20);

	TDetM20 = _M10 * (_M21 * _M33 - _M23 * _M31) -
		_M11 * (_M20 * _M33 - _M23 * _M30) +
		_M13 * (_M20 * _M31 - _M21 * _M30);

	TDetM21 = _M00 * (_M21 * _M33 - _M23 * _M31) -
		_M01 * (_M20 * _M33 - _M23 * _M30) +
		_M03 * (_M20 * _M31 - _M21 * _M30);

	TDetM22 = _M00 * (_M11 * _M33 - _M13 * _M31) -
		_M01 * (_M10 * _M33 - _M13 * _M30) +
		_M03 * (_M10 * _M31 - _M11 * _M30);

	TDetM23 = _M00 * (_M11 * _M23 - _M13 * _M21) -
		_M01 * (_M10 * _M23 - _M13 * _M20) +
		_M03 * (_M10 * _M21 - _M11 * _M20);

	TDetM30 = _M10 * (_M21 * _M32 - _M22 * _M31) -
		_M11 * (_M20 * _M32 - _M22 * _M30) +
		_M12 * (_M20 * _M31 - _M21 * _M30);

	TDetM31 = _M00 * (_M21 * _M32 - _M22 * _M31) -
		_M01 * (_M20 * _M32 - _M22 * _M30) +
		_M02 * (_M20 * _M31 - _M21 * _M30);

	TDetM32 = _M00 * (_M11 * _M32 - _M12 * _M31) -
		_M01 * (_M10 * _M32 - _M12 * _M30) +
		_M02 * (_M10 * _M31 - _M11 * _M30);

	TDetM33 = _M00 * (_M11 * _M22 - _M12 * _M21) -
		_M01 * (_M10 * _M22 - _M12 * _M20) +
		_M02 * (_M10 * _M21 - _M11 * _M20);

	_M00 = TDetM00 * aDet; _M01 = -TDetM01 * aDet; _M02 = TDetM02 * aDet; _M03 = -TDetM03 * aDet;
	_M10 = -TDetM10 * aDet; _M11 = TDetM11 * aDet; _M12 = -TDetM12 * aDet; _M13 = TDetM13 * aDet;
	_M20 = TDetM20 * aDet; _M21 = -TDetM21 * aDet; _M22 = TDetM22 * aDet; _M23 = -TDetM23 * aDet;
	_M30 = -TDetM30 * aDet; _M31 = TDetM31 * aDet; _M32 = -TDetM32 * aDet; _M33 = TDetM33 * aDet;

	return *this;
}

// 產生轉置矩陣
CMatrix4FLOAT& CMatrix4FLOAT::Transpose(void) {
	*this = CMatrix4FLOAT(_M00, _M01, _M02, _M03,
		_M10, _M11, _M12, _M13,
		_M20, _M21, _M22, _M23,
		_M30, _M31, _M32, _M33);
	return *this;
}

// 覆載 + 法運算子 Mo = M + Mi
CMatrix4FLOAT CMatrix4FLOAT::operator+(const CMatrix4FLOAT & Matrix) const {
	return Add(Matrix);
}

// 覆載 - 法運算子 Mo = M - Mi
CMatrix4FLOAT CMatrix4FLOAT::operator-(const CMatrix4FLOAT & Matrix) const {
	return Sub(Matrix);
}

// 覆載 * 法對實數的運算子	Mo = M * Value
CMatrix4FLOAT CMatrix4FLOAT::operator*(const float Value) const {
	return Mul(Value);
}

// 覆載 * 法對矩陣的運算子 Mo = M * Mi
CMatrix4FLOAT CMatrix4FLOAT::operator*(const CMatrix4FLOAT & Matrix) const {
	return Mul(Matrix);
}

// 矩陣與向量的乘法 [m_x, m_y, m_z, 1] = M * [m_x, m_y, m_z, 1]
CVector3FLOAT CMatrix4FLOAT::operator*(const CVectorReference3FLOAT & Vector) const {
	return Mul(Vector);
}

// 矩陣與向量的乘法 V = M * Vi
CVector4FLOAT CMatrix4FLOAT::operator*(const CVectorReference4FLOAT & Vector) const {
	return Mul(Vector);
}

// 覆載 / 法對實數的運算子 Mo = M / Value
CMatrix4FLOAT CMatrix4FLOAT::operator/(const float Value) const {
	return Div(Value);
}

// 使兩矩陣相等 M = Mi
CMatrix4FLOAT& CMatrix4FLOAT::operator=(const CMatrix4FLOAT & Matrix) {
	Assign(Matrix);
	return *this;
}

// 覆載 += 法運算子 M = M + Mi
CMatrix4FLOAT& CMatrix4FLOAT::operator+=(const CMatrix4FLOAT & Matrix) {
	AddAssign(Matrix);
	return *this;
}

// 覆載 -= 法運算子 M = M - Mi
CMatrix4FLOAT& CMatrix4FLOAT::operator-=(const CMatrix4FLOAT & Matrix) {
	SubAssign(Matrix);
	return *this;
}

// 覆載 *= 法運算子 M = M * Value
CMatrix4FLOAT& CMatrix4FLOAT::operator*=(const float Value) {
	MulAssign(Value);
	return *this;
}

// 覆載 *= 法對矩陣的運算子 M = M * Mi
CMatrix4FLOAT& CMatrix4FLOAT::operator*=(const CMatrix4FLOAT & Matrix) {
	MulAssign(Matrix);
	return *this;
}

// 覆載 /= 法對實數的運算子 M = M / Value
CMatrix4FLOAT& CMatrix4FLOAT::operator/=(const float Value) {
	DivAssign(Value);
	return *this;
}

// 使兩矩陣相加 Mo = M + Mi
CMatrix4FLOAT CMatrix4FLOAT::Add(const CMatrix4FLOAT & Matrix) const {
	return CMatrix4FLOAT(	_M00 + Matrix._M00, _M10 + Matrix._M10, _M20 + Matrix._M20, _M30 + Matrix._M30,
							_M01 + Matrix._M01, _M11 + Matrix._M11, _M21 + Matrix._M21, _M31 + Matrix._M31,
							_M02 + Matrix._M02, _M12 + Matrix._M12, _M22 + Matrix._M22, _M32 + Matrix._M32,
							_M03 + Matrix._M03, _M13 + Matrix._M13, _M23 + Matrix._M23, _M33 + Matrix._M33);
}

// 使兩矩陣相減 Mo = M - Mi
CMatrix4FLOAT CMatrix4FLOAT::Sub(const CMatrix4FLOAT & Matrix) const {
	return CMatrix4FLOAT(	_M00 - Matrix._M00, _M10 - Matrix._M10, _M20 - Matrix._M20, _M30 - Matrix._M30,
							_M01 - Matrix._M01, _M11 - Matrix._M11, _M21 - Matrix._M21, _M31 - Matrix._M31,
							_M02 - Matrix._M02, _M12 - Matrix._M12, _M22 - Matrix._M22, _M32 - Matrix._M32,
							_M03 - Matrix._M03, _M13 - Matrix._M13, _M23 - Matrix._M23, _M33 - Matrix._M33);

}

//  矩陣與實數做乘法 Mo = M * Value
CMatrix4FLOAT CMatrix4FLOAT::Mul(const float Value) const {
	return CMatrix4FLOAT(	_M00 * Value, _M10 * Value, _M20 * Value, _M30 * Value,
							_M01 * Value, _M11 * Value, _M21 * Value, _M31 * Value,
							_M02 * Value, _M12 * Value, _M22 * Value, _M32 * Value,
							_M03 * Value, _M13 * Value, _M23 * Value, _M33 * Value);
}

// 矩陣的乘法 Mo = M * Mi
CMatrix4FLOAT CMatrix4FLOAT::Mul(const CMatrix4FLOAT & Matrix) const {
	return CMatrix4FLOAT(	_M00 * Matrix._M00 + _M01 * Matrix._M10 + _M02 * Matrix._M20 + _M03 * Matrix._M30, // M00
							_M10 * Matrix._M00 + _M11 * Matrix._M10 + _M12 * Matrix._M20 + _M13 * Matrix._M30, // M10
							_M20 * Matrix._M00 + _M21 * Matrix._M10 + _M22 * Matrix._M20 + _M23 * Matrix._M30, // M20
							_M30 * Matrix._M00 + _M31 * Matrix._M10 + _M32 * Matrix._M20 + _M33 * Matrix._M30, // M30
							_M00 * Matrix._M01 + _M01 * Matrix._M11 + _M02 * Matrix._M21 + _M03 * Matrix._M31, // M01
							_M10 * Matrix._M01 + _M11 * Matrix._M11 + _M12 * Matrix._M21 + _M13 * Matrix._M31, // M11
							_M20 * Matrix._M01 + _M21 * Matrix._M11 + _M22 * Matrix._M21 + _M23 * Matrix._M31, // M21
							_M30 * Matrix._M01 + _M31 * Matrix._M11 + _M32 * Matrix._M21 + _M33 * Matrix._M31, // M31
							_M00 * Matrix._M02 + _M01 * Matrix._M12 + _M02 * Matrix._M22 + _M03 * Matrix._M32, // M02
							_M10 * Matrix._M02 + _M11 * Matrix._M12 + _M12 * Matrix._M22 + _M13 * Matrix._M32, // M12
							_M20 * Matrix._M02 + _M21 * Matrix._M12 + _M22 * Matrix._M22 + _M23 * Matrix._M32, // M22
							_M30 * Matrix._M02 + _M31 * Matrix._M12 + _M32 * Matrix._M22 + _M33 * Matrix._M32, // M32
							_M00 * Matrix._M03 + _M01 * Matrix._M13 + _M02 * Matrix._M23 + _M03 * Matrix._M33, // M03
							_M10 * Matrix._M03 + _M11 * Matrix._M13 + _M12 * Matrix._M23 + _M13 * Matrix._M33, // M13
							_M20 * Matrix._M03 + _M21 * Matrix._M13 + _M22 * Matrix._M23 + _M23 * Matrix._M33, // M23
							_M30 * Matrix._M03 + _M31 * Matrix._M13 + _M32 * Matrix._M23 + _M33 * Matrix._M33);// M33
}

// 矩陣與向量的乘法 [m_x, m_y, m_z, 1] = M * [m_x, m_y, m_z, 1]
CVector3FLOAT CMatrix4FLOAT::Mul(const CVectorReference3FLOAT & Vector) const {
	float w = _M30 * Vector.m_x + _M31 * Vector.m_y + _M32 * Vector.m_z + _M33;
	if (w == 0) {
		fprintf(stderr, "\n error:Matrix4f Mul Vec3f w == 0!.");
		w = 1;
	}
	float aw = 1.0f / w;

	return CVector3FLOAT(	(_M00 * Vector.m_x + _M01 * Vector.m_y + _M02 * Vector.m_z + _M03) * aw,
							(_M10 * Vector.m_x + _M11 * Vector.m_y + _M12 * Vector.m_z + _M13) * aw,
							(_M20 * Vector.m_x + _M21 * Vector.m_y + _M22 * Vector.m_z + _M23) * aw);
}

// 矩陣與向量的乘法 V = M * Vi
CVector4FLOAT CMatrix4FLOAT::Mul(const CVectorReference4FLOAT & Vector) const {
	return CVector4FLOAT(	_M00 * Vector.m_x + _M01 * Vector.m_y + _M02 * Vector.m_z + _M03 * Vector.m_w,
							_M10 * Vector.m_x + _M11 * Vector.m_y + _M12 * Vector.m_z + _M13 * Vector.m_w,
							_M20 * Vector.m_x + _M21 * Vector.m_y + _M22 * Vector.m_z + _M23 * Vector.m_w,
							_M30 * Vector.m_x + _M31 * Vector.m_y + _M32 * Vector.m_z + _M33 * Vector.m_w);

}

//  矩陣與實數做除法 Mo = M / Value
CMatrix4FLOAT CMatrix4FLOAT::Div(const float Value) const {
	if (Value == 0) {
		fprintf(stderr, "\n error:Matrix4f Div Value == 0!.");
		return *this;
	}
	return CMatrix4FLOAT(	_M00 / Value, _M10 / Value, _M20 / Value, _M30 / Value,
							_M01 / Value, _M11 / Value, _M21 / Value, _M31 / Value,
							_M02 / Value, _M12 / Value, _M22 / Value, _M32 / Value,
							_M03 / Value, _M13 / Value, _M23 / Value, _M33 / Value);
}

// 使兩矩陣相等 M = Mi
void CMatrix4FLOAT::Assign(const CMatrix4FLOAT & Matrix) {
	_M00 = Matrix._M00; _M01 = Matrix._M01; _M02 = Matrix._M02; _M03 = Matrix._M03;
	_M10 = Matrix._M10; _M11 = Matrix._M11; _M12 = Matrix._M12; _M13 = Matrix._M13;
	_M20 = Matrix._M20; _M21 = Matrix._M21; _M22 = Matrix._M22; _M23 = Matrix._M23;
	_M30 = Matrix._M30; _M31 = Matrix._M31; _M32 = Matrix._M32; _M33 = Matrix._M33;
}

// 計算兩矩陣相加並且保存數值 M = M + Mi
void CMatrix4FLOAT::AddAssign(const CMatrix4FLOAT & Matrix) {
	_M00 += Matrix._M00; _M01 += Matrix._M01; _M02 += Matrix._M02; _M03 += Matrix._M03;
	_M10 += Matrix._M10; _M11 += Matrix._M11; _M12 += Matrix._M12; _M13 += Matrix._M13;
	_M20 += Matrix._M20; _M21 += Matrix._M21; _M22 += Matrix._M22; _M23 += Matrix._M23;
	_M30 += Matrix._M30; _M31 += Matrix._M31; _M32 += Matrix._M32; _M33 += Matrix._M33;
}

// 計算兩矩陣相減並且保存數值 M = M - Mi
void CMatrix4FLOAT::SubAssign(const CMatrix4FLOAT & Matrix) {
	_M00 -= Matrix._M00; _M01 -= Matrix._M01; _M02 -= Matrix._M02; _M03 -= Matrix._M03;
	_M10 -= Matrix._M10; _M11 -= Matrix._M11; _M12 -= Matrix._M12; _M13 -= Matrix._M13;
	_M20 -= Matrix._M20; _M21 -= Matrix._M21; _M22 -= Matrix._M22; _M23 -= Matrix._M23;
	_M30 -= Matrix._M30; _M31 -= Matrix._M31; _M32 -= Matrix._M32; _M33 -= Matrix._M33;
}

// 計算兩矩陣相乘並且保存數值 M = M * Mi
void CMatrix4FLOAT::MulAssign(const CMatrix4FLOAT & Matrix) {
	*this = (*this).Mul(Matrix);
}

//  矩陣與實數做乘法並且保存數值 M = M * Value
void CMatrix4FLOAT::MulAssign(const float Value) {
	_M00 *= Value;	_M10 *= Value;	_M20 *= Value;	_M30 *= Value;
	_M01 *= Value;	_M11 *= Value;	_M21 *= Value;	_M31 *= Value;
	_M02 *= Value;	_M12 *= Value;	_M22 *= Value;	_M32 *= Value;
	_M03 *= Value;	_M13 *= Value;	_M23 *= Value;	_M33 *= Value;
}

//  矩陣與實數做除法並且保存數值 M = M / Value
void CMatrix4FLOAT::DivAssign(const float Value) {
	if (Value == 0) {
		fprintf(stderr, "\n error:Matrix4f DivAssign Value == 0!.");
		return;
	}
	_M00 /= Value;	_M10 /= Value;	_M20 /= Value;	_M30 /= Value;
	_M01 /= Value;	_M11 /= Value;	_M21 /= Value;	_M31 /= Value;
	_M02 /= Value;	_M12 /= Value;	_M22 /= Value;	_M32 /= Value;
	_M03 /= Value;	_M13 /= Value;	_M23 /= Value;	_M33 /= Value;
}
