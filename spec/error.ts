export class ResponseCodedError extends Error {
  constructor(public error_code: number, message?: string) {
    super(message)
    Object.setPrototypeOf(this, ResponseCodedError.prototype);
  }
  __ResponseCodedError = true
  static isinstance(instance: unknown) {
    return typeof instance === 'object'
      && instance !== null
      && '__ResponseCodedError' in instance
      && instance.__ResponseCodedError === true
  }
}

/**
 * The method is not supported
 */
export class UnsupportedMethodError extends ResponseCodedError {
  constructor() {
    super(405, 'Unsupported method')
    Object.setPrototypeOf(this, UnsupportedMethodError.prototype)
  }
  __UnsupportedMethodError = true
  static isinstance(instance: unknown) {
    return typeof instance === 'object'
      && instance !== null
      && '__UnsupportedMethodError' in instance
      && instance.__UnsupportedMethodError === true
  }
}

/**
 * The resource is not found
 */
export class NotFoundError extends ResponseCodedError {
  constructor() {
    super(404, 'Not Found')
    Object.setPrototypeOf(this, NotFoundError.prototype)
  }
  __NotFoundError = true
  static isinstance(instance: unknown) {
    return typeof instance === 'object'
      && instance !== null
      && '__NotFoundError' in instance
      && instance.__NotFoundError === true
  }
}

/**
 * The resource is not found
 */
export class UnauthorizedError extends ResponseCodedError {
  constructor() {
    super(401, 'Unauthorized')
    Object.setPrototypeOf(this, UnauthorizedError.prototype)
  }
  __UnauthorizedError = true
  static isinstance(instance: unknown) {
    return typeof instance === 'object'
      && instance !== null
      && '__UnauthorizedError' in instance
      && instance.__UnauthorizedError === true
  }
}

/**
 * This timeout error is used to ensure we don't wait too long for dependencies
 *  fortunately, even if it occurs the job will requeue still making progress.
 */
export class TimeoutError extends ResponseCodedError {
  constructor() {
    super(504, 'Timeout reached')
    Object.setPrototypeOf(this, TimeoutError.prototype)
  }
  __TimeoutError = true
  static isinstance(instance: unknown) {
    return typeof instance === 'object'
      && instance !== null
      && '__TimeoutError' in instance
      && instance.__TimeoutError === true
  }
}

/**
 * This error occurs when the input node is not populated yet
 */
export class UnboundError extends ResponseCodedError {
  constructor() {
    super(422, 'Refusing to submit unbound variable')
    Object.setPrototypeOf(this, UnboundError.prototype)
  }
  __UnboundError = true
  static isinstance(instance: unknown) {
    return typeof instance === 'object'
      && instance !== null
      && '__UnboundError' in instance
      && instance.__UnboundError === true
  }
}
