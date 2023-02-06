export class ResponseCodedError extends Error {
  constructor(public error_code: number, message?: string) {
    super(message)
    Object.setPrototypeOf(this, ResponseCodedError.prototype);
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
}

/**
 * The resource is not found
 */
export class NotFoundError extends ResponseCodedError {
  constructor() {
    super(404, 'Not Found')
    Object.setPrototypeOf(this, NotFoundError.prototype)
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
}

/**
 * This error occurs when the input node is not populated yet
 */
export class UnboundError extends ResponseCodedError {
  constructor() {
    super(422, 'Refusing to submit unbound variable')
    Object.setPrototypeOf(this, UnboundError.prototype)
  }
}
